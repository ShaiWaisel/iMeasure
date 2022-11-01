%close all; clear all;clc; imtool close all;
VERBOSE = 1;
NOVERBOSE = 0;
%read frames for all cameras
platformID= "BaseP001";
sensor1ID = 'f1381564';
sensor2ID = 'f1371657';
sensor3ID = 'f1230436';
sensor4ID = 'f0000000';

connPlatform = sqlite('C:\Projects\iMeasure\FromErez\Database\!iMeasure.db'); 
platforms = sqlread(connPlatform,'Platforms');

sqlQuery = sprintf('SELECT * FROM Platforms WHERE PlatformID = "%s"', platformID);
newdata =['UPDATE Platforms' ,sprintf(' SET Sensor1_DeviceID = "%s", Sensor2_DeviceID = "%s", Sensor3_DeviceID = "%s", Sensor4_DeviceID = "%s"',...
    sensor1ID,sensor2ID,sensor3ID, sensor4ID), sprintf(' WHERE PlatformID = "%s"', platformID)];
exec(connPlatform, newdata);

cameras = [sensor1ID; sensor2ID; sensor3ID];

conn = sqlite('C:\Projects\iMeasure\FromErez\Database\SamCore11.db'); 
for camera=1:size(cameras,1)
    rawData{camera} = sqlread(conn,sprintf('%s_Frames',cameras(camera,:)));
    wall{camera}= [];
    transWall{camera}= [];
    Uval{camera} = [];
    Vval{camera} = [];
    moving{camera} = [];
    fixed{camera} = [];
    transMats{camera} = [];
end
close(conn);
% set global calibration
cameraPosition = [0.7, 0.6, 0.5];
lookAtPoint = [1.5, 1.5, 0.5];
squareSize = 0.1;
boardSize = 1.5;

%create theoretical grid (centers of squares)
[worldGrid.rGrid, worldGrid.gGrid, worldGrid.bGrid, worldGrid.yGrid] = ...
    WorldGrid(squareSize, boardSize, cameraPosition, lookAtPoint, NOVERBOSE);

FRAMES = min(size(rawData{1}.ID,1),8);
U=[];
V=[];
tform = affine3d;
loop = 0;
for frame=1:FRAMES
    for camera=1:size(cameras,1)
        [original, u, v, img] = GrabFrameData(rawData{camera}, frame);
        Uval{camera} = [Uval{camera}; u];
        Vval{camera} = [Vval{camera}; v];
        if (frame == 1)
            wall{camera} = original;
        else
            wall{camera} = pccat([wall{camera}, original]);
        end
    end
end
if (VERBOSE)
    hSegmentation = figure; axis equal; hold on; 
    hAxes = gca;
end;
for camera=1:size(cameras,1)
    [~, ~, ~, image{camera}] = GrabFrameData(rawData{camera},1);
    for colorCode=1:4
        if (size(transMats{camera},1) >= 0)
            [code, ROI, mask] = DetectColors(image{camera},colorCode, NOVERBOSE);
            if (code == 1)
                % found colorCode area
                squares =  AnalyzeBoard3(image{camera}, ROI, mask, colorCode, NOVERBOSE);
                for i=1:size(squares,1)
                    for j=1:size(squares,2)
                        if (squares(i,j).value > 0)
                            center2D = mean(squares(i,j).quad);
                            featuresIdxs = Features(Uval{camera}, Vval{camera}, image{camera}, center2D);
                            if (size(featuresIdxs,1) > 0)
                                center3D = wall{camera}.Location(featuresIdxs,:);
                                moving{camera} = [moving{camera}; center3D];
                                switch colorCode
                                    case 1
                                        fixed{camera} = [fixed{camera}; reshape(worldGrid.rGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                    case 2 
                                        fixed{camera} = [fixed{camera}; reshape(worldGrid.gGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                    case 3 
                                        fixed{camera} = [fixed{camera}; reshape(worldGrid.bGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                    case 4 
                                        fixed{camera} = [fixed{camera}; reshape(worldGrid.yGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    m=moving{camera};
    f=fixed{camera};
    if (size(m,1)>2)
        [mat, r, e] = absor(m', f')
        f1=[f, ones(size(f,1),1)];
        m1=[m, ones(size(m,1),1)];
        transMats{camera} = mat.M';
        z=m1*mat.M';
        if (NOVERBOSE)
            hold on; scatter3(hAxes, f(:,1), f(:,2), f(:,3),'ob'); scatter3(hAxes, m(:,1), m(:,2), m(:,3),'+r');
            for i=1:size(f,1)
                text(hAxes, f(i,1),f(i,2),f(i,3)+50,num2str(i));
                    text(hAxes, m(i,1),m(i,2),m(i,3)+50,num2str(i));
            end
            scatter3(hAxes, z(:,1), z(:,2), z(:,3),'*k')
            for i=1:size(f,1)
                text(hAxes, f(i,1),f(i,2),f(i,3)+50,num2str(i),'Color','blue');
                text(hAxes, m(i,1),m(i,2),m(i,3)+50,num2str(i),'Color','red');
                text(hAxes, z(i,1),z(i,2),z(i,3)+50,num2str(i),'Color','black');
            end
        end
    end

end
%[planes, modified, planars] = SetPlanarAreas(rightWall, 0.005, 300000);
%rightWall = pcdownsample(modified, 'gridAverage', 0.0005);
for camera=1:size(cameras,1)
    sqSensorName = sprintf('Sensor%d_tForm', camera);
    mat = transMats{camera}; 
    sqMat='';
    for i=1:4
        for j=1:4
            sqMat = strcat(sqMat, sprintf('%f; ',mat(i,j)));
        end
    end
    sqlOp =['UPDATE Platforms' , ...
        sprintf(' SET %s = "%s" ', sqSensorName, sqMat), ...
        sprintf(' WHERE PlatformID = "%s"', platformID)];
    exec(connPlatform, sqlOp);
end
cloud4 = pointCloud([-1,1 -1;1 1 -1; -1 1 1; 1 1 1]);
close(connPlatform);
for camera=1:size(cameras,1)
    mat = transMats{camera}; 
    tform = affinetform3d(mat');

    transWall{camera} = pctransform(wall{camera}, tform);
    %transWall{camera} = pctransform(cloud4, tform);
end
close(conn);
player = pcplayer([-1.5 1.5], [-1.5 1.5], [0 1]);
walls = pccat([transWall{1} transWall{2} transWall{3}]);
view(player, walls);
return;
figure;
axis equal;
hold on;
for camera=1:size(cameras,1)
    hold on;
    pcshow(transWall{camera});
end
return;


% ===================================================================================================

function [ptc, fu, fv, image] = GrabFrameData(data, frame)
    rgbW = 1280;
    rgbH = 720;
    depthW = 640;
    depthH = 480;
    rgb = data(frame,19);
    rgb = cell2mat(rgb.FrameColor);
    rgb = reshape(rgb,3,[]);
    r = reshape(rgb(1,:), rgbW, rgbH)'; 
    g = reshape(rgb(2,:), rgbW, rgbH)'; 
    b = reshape(rgb(3,:), rgbW, rgbH)';
    image = cat(3, r, g, b);
    %figure, imshow(image);
    x = (cell2mat(data(frame, 20).FrameX));
    y = (cell2mat(data(frame, 21).FrameY));
    z = (cell2mat(data(frame, 22).FrameZ));
    u = (cell2mat(data(frame, 23).FrameU));
    v = (cell2mat(data(frame, 24).FrameV));
    fx=typecast(x,'single');
    fy=typecast(y,'single');
    fz=typecast(z,'single');
    fu=typecast(u,'single');
    fv=typecast(v,'single');
    idx = (fu<0.1) | (fu>0.9) | (fv<0.1) | (fv>0.9);
    fx(idx)=[];
    fy(idx)=[];
    fz(idx) = [];
    fu(idx) = [];
    fv(idx) = [];
    iu = int32(floor(size(image,2)*fu)+1);
    iv = int32(floor(size(image,1)*fv)+1);
    linIdx = sub2ind(size(image), iv, iu);
    ptc = pointCloud([fx, fz, -fy]);
    ptc.Color=[r(linIdx), g(linIdx), b(linIdx)];
end

%-------------------------------------------------------------------------------------------------

function [moving, fixed] = DetectCheckerboardRight(image, wallCloud, u, v, centerIdx, fixedAngle, originPosition, cameraPosition)
[blackPoints, colorPoints, centerPoint] = DetectSquareCorners(image);
featuresIdxs = Features(u, v, image, centerPoint);
centerWall = wallCloud.Location(featuresIdxs,:);
featuresIdxs = Features(u, v, image, colorPoints);
moving = wallCloud.Location(featuresIdxs,:);
mask = zeros(11,11);
for idx=1:size(moving,1)
    testLoc = moving(idx,:);
    distXY = norm(testLoc(1:2) - centerWall(1:2));
    distZ = centerWall(3) - testLoc(3);
    s = centerIdx + sign(testLoc(1) - centerWall(1)) * round(distXY*10); 
    mask(centerIdx + round(distZ*10), centerIdx + sign(testLoc(1) - centerWall(1)) * round(distXY*10)) = 1;
end
fixed = [];
mask(centerIdx,centerIdx) = 2;
for i=1:11
    for j=1:11
        if (mask(i,j) > 0)
            X = (j-1) * 0.1 * cosd(fixedAngle);
            Y = norm(cameraPosition - originPosition)' - (j-1) * 0.1 * sind(fixedAngle);
            Z = cameraPosition(3) - (i-1) * 0.1;
            fixed = [fixed; X, Y, Z];

        end
    end
end
fixed = pointCloud(fixed);
moving = pointCloud(moving);
end

%-------------------------------------------------------------------------------------------------

function [moving, fixed] = DetectCheckerboardLeft(image, wallCloud, u, v, centerIdx, fixedAngle, originPosition, cameraPosition)
[blackPoints, colorPoints, centerPoint] = DetectSquareCorners(image);
featuresIdxs = Features(u, v, image, centerPoint);
centerWall = wallCloud.Location(featuresIdxs,:);
featuresIdxs = Features(u, v, image, colorPoints);
moving = wallCloud.Location(featuresIdxs,:);
mask = zeros(11,11);
for idx=1:size(moving,1)
    testLoc = moving(idx,:);
    distXY = norm(testLoc(1:2) - centerWall(1:2));
    distZ = centerWall(3) - testLoc(3);
    s = centerIdx + sign(testLoc(1) - centerWall(1)) * round(distXY*10); 
    mask(centerIdx + round(distZ*10), centerIdx + sign(testLoc(1) - centerWall(1)) * round(distXY*10)) = 1;
end
fixed = [];
mask(centerIdx,centerIdx) = 2;
for i=1:11
    for j=1:11
        if (mask(i,j) > 0)
            X = (1-j) * 0.1 * cosd(fixedAngle);
            Y = norm(cameraPosition - originPosition)' - (j-1) * 0.1 * sind(fixedAngle);
            Z = cameraPosition(3) - (i-1) * 0.1;
            fixed = [fixed; X, Y, Z];

        end
    end
end
fixed = pointCloud(fixed);
moving = pointCloud(moving);
end

%-------------------------------------------------------------------------------------------------

function [planes, modified, planars, boundsList, edges, corners] = SetPlanarAreas(ptc, maxDistance, planeInPixels)
    referenceVector = [0,1,0];
    maxAngularDistance = 80;
    modified = 0;
    boundsList=[];
    intersectionsList = [];
    planes2go = 1;      
    planes = [];
    planars = pointCloud([0,0,0]);
    while (planes2go > 0)
        [plane,inlierIndices,outlierIndices] = pcfitplane(ptc, maxDistance,referenceVector,maxAngularDistance);
        if (size(inlierIndices, 1) > planeInPixels)
            planar = select(ptc,inlierIndices);
            distFromPlane = plane.Parameters(1)*planar.Location(:,1) + ...
                plane.Parameters(2)*planar.Location(:,2) + ...
                plane.Parameters(3)*planar.Location(:,3) + ...
                plane.Parameters(4);
            newLocs = planar.Location - plane.Normal .* distFromPlane;   
            newPlanar = pointCloud(newLocs);
            planars = pccat([planars, newLocs]);
            outliers = select(ptc,outlierIndices);
            outliers.Color = ptc.Color(outlierIndices,:);
            newCloud = pointCloud(newLocs);
            planes = [planes; [plane.Parameters, newCloud.Location(1,:)] ];
            newCloud.Color = ptc.Color(inlierIndices, :);
            if (modified == 0)
                modified = newCloud;
            else
                modified = pccat([modified; newCloud]);
            end
            ptc = outliers;
        else
            planes2go = 0;
        end
    end
    modified = pccat([modified; ptc]);
    intersects = [];
    for i=1:size(planes,1)
        for j=i+1:size(planes,1)
            [P, N, check] = plane_intersect(planes(i,1:3), planes(i,5:7), planes(j,1:3), planes(j,5:7));
            intersects = [intersects; P, N];
        end
    end
    corners = [];
    edges = [];
    span = 3;
    for i=1:size(intersects,1)
        P0 = intersects(i,1:3)-intersects(i,4:6) * span;
        P1 = intersects(i,1:3)+intersects(i,4:6) * span;
        edges = [edges; P0; P1;[NaN, NaN, NaN]];
        for j=i+1:size(intersects,1)
            Q0 = intersects(j,1:3)-intersects(j,4:6) * span;
            Q1 = intersects(j,1:3)+intersects(j,4:6) * span;
            [D,Xcp,Ycp,Zcp,Xcq,Ycq,Zcq,Dmin,imin,jmin]= ll_dist3d(P0, P1, Q0, Q1);
            if (D < 0.01) && (Xcp > -span) && (Xcp < span) && (Ycp > -span) && (Ycp < span) && (Zcp > -span) && (Zcp < span)
                corners = [corners; Xcp, Ycp, Zcp];
            end
        end
    end
end

%-------------------------------------------------------------------------------------------------

function [edges, corners] = transformEdgesCorners(e, c, tform)
e = [e, ones(sizeof(e,1))];
edges = e * tform.T;
c = [c, ones(sizeof(c,1))];
corners = c * tform.T;
end

%-------------------------------------------------------------------------------------------------

function [bounds] = PlaneBounds(pts, pe)
    idx = randi([1, size(pts,1)],1000,1);
    b = boundary(pts(idx,1), pts(idx,3));
    ib = idx(b);
    bpy = -(pe(1)*pts(ib,1) + pe(3)*pts(ib,3) + pe(4))/pe(2);
    bounds = [pts(ib,1), bpy, pts(ib,3)];
end

%-------------------------------------------------------------------------------------------------

function [unified] = UnifyExtended(baseCloud, newCloud)
    outliers = ((newCloud.Location(:,1) < baseCloud.XLimits(1)) | ...
                (newCloud.Location(:,1) > baseCloud.XLimits(2)) | ...
                (newCloud.Location(:,2) < baseCloud.YLimits(1)) | ...
                (newCloud.Location(:,2) > baseCloud.YLimits(2)) | ...
                (newCloud.Location(:,3) < baseCloud.ZLimits(1)) | ...
                (newCloud.Location(:,3) > baseCloud.ZLimits(2))) ;
    newCloud = select(newCloud, outliers);
    unified = pccat([baseCloud, newCloud]);
    
end

%-------------------------------------------------------------------------------------------------

function [featuresIndices] = Features(cloudU, cloudV, image, featuersPos)

    fuv = [featuersPos(:,1) / size(image,2), featuersPos(:,2) / size(image,1)];
    featuresIndices = [];
    for idx=1:size(fuv,1)
        if (fuv(idx,1) > 0.2) & (fuv(idx,1) < 0.8) & (fuv(idx,2) > 0.2) & (fuv(idx,2) < 0.8)

            dists = (fuv(idx, 1) - cloudU).*(fuv(idx,1) -cloudU) + (fuv(idx, 2) - cloudV).*(fuv(idx,2) -cloudV);
        
            [d,idx] = min(dists);
            featuresIndices = [featuresIndices, idx];
        end
    end
end

%-------------------------------------------------------------------------------------------------

function [blackCorners, colorCorners, centerPoint] = DetectSquareCorners(image)

[xnd,map] = rgb2ind(image,3,'nodither');

img = uint8(rgb2gray(ind2rgb(xnd, map))*255);
idxs = sort(find(imhist(img)>0)-1);
grays = (img == idxs(2));
checkerBoard = img;
checkerBoard(grays) = idxs(1);
[colorCorners, boardSize] = detectCheckerboardPoints(checkerBoard);
checkerBoard(grays) = idxs(3);
[blackCorners, boardSize] = detectCheckerboardPoints(checkerBoard);
centerPoint=[];
if (boardSize == [4, 4])
    ctr = [mean(blackCorners)];
    [v,idx]=min((blackCorners(:,1) - ctr(1)).^2 + (blackCorners(:,2) - ctr(2)).^2);
    centerPoint = blackCorners(idx,:);
end

figure, imshow(image);
hold on;
scatter(colorCorners(:,1), colorCorners(:,2),'o','filled');
for i=1:size(colorCorners,1)
    b = num2str(i); c=cellstr(b);
    text(colorCorners(i,1)+10, colorCorners(i,2)+10, c,'FontSize', 20, 'Color','magenta');
end
    
end







