close all; clear all;
conn = sqlite('C:\Projects\iMeasure\FromErez\Database\H001_Calib_22_11_22.db'); 
leftData = sqlread(conn,'f1381564_Frames');
rightData = sqlread(conn,'f1371657_Frames');
backData = sqlread(conn,'f1230436_Frames');
close(conn);
cameraPosition = [700, 600, 500];
lookAtPoint = [1500, 1500, 500];
squareSize = 100;
%worldGrid = WorldGrid(squareSize, cameraPosition, lookAtPoint);

FRAMES = min(min(size(rightData,1),size(leftData,1)),8);
cameraPosition = [0.7, 0.7, 0.5];
originPosition = [0.0, 0.0, 0.5];
rightAngle = atan2d(cameraPosition(2) - originPosition(2), cameraPosition(1) - originPosition(1));
s = [0:0.1:1]';
X = s * cosd(rightAngle);
Y = norm(cameraPosition - originPosition)' - s * sind(rightAngle);
Z = [1:-0.1:0]' - cameraPosition(3);
rightPoints = [X, Y, Z];
leftAngle = 90.0 - rightAngle;
X = -s * sind(leftAngle);
Y = norm(cameraPosition - originPosition) - s * cosd(leftAngle);
leftPoints = [X, Y, Z];
centerIdx = 6;

gridSize = 0.01;
nonuniformGridSample = 6;
h2=[];
rightWall = pointCloud([0,0,0]);
fullFixed = rightWall;
scannedModel = rightWall;
accumTform = affine3d;
transMats = affine3d;
model = [];
tform = affine3d;
loop = 0;
[rightWall,rightU, rightV, rightImage] = GrabFrameData(rightData,1); 
[leftWall,leftU, leftV, leftImage] = GrabFrameData(leftData,1); 
[backWall,backU, backV, backImage] = GrabFrameData(backData,1); 
for colorCode=1:4
    [code, croppedImage, croppedMask] = DetectColors(rightImage,colorCode);
    if (code == 1)
        squares =  AnalyzeBoard3(croppedImage, croppedMask, colorCode);

    end
end
return;
for frame=2:1:FRAMES
    [original, u, v, rightImage] = GrabFrameData(rightData, frame);
    rightU = [rightU; u];
    rightV = [rightV; v];
    rightWall = pccat([rightWall, original]);
    [original, u, v, leftImage] = GrabFrameData(leftData, frame);
    leftU = [leftU; u];
    leftV = [leftV; v];
    leftWall = pccat([leftWall, original]);
    [original, u, v, backImage] = GrabFrameData(backData, frame);
    backU = [backU; u];
    backV = [backV; v];
    leftWall = pccat([leftWall, original]);
end
[movingRight, fixedRight] = DetectCheckerboardRight(rightImage, rightWall, rightU, rightV, centerIdx, rightAngle, originPosition, cameraPosition)

%tform = pcregistericp(moving, fixed, 'Tolerance', [0.1, 0.05], 'InitialTransform', rigid3d(rotz(30), [0,0,0]));
tformRight = pcregisterndt(movingRight, fixedRight, 1.414);
rightWall = pctransform(rightWall, tformRight);

[movingLeft, fixedLeft] = DetectCheckerboardLeft(leftImage, leftWall, leftU, leftV, centerIdx, leftAngle, originPosition, cameraPosition)
tformLeft = pcregisterndt(movingLeft, fixedLeft, 1.414);
leftWall = pctransform(leftWall, tformLeft);


% figure, imshow(rightImage);
% hold on;
% scatter(points(:,1), points(:,2),'yo','filled');
% featuresIndices = Features()

[planes, modified, planars] = SetPlanarAreas(rightWall, 0.005, 300000);
%rightWall = pcdownsample(modified, 'gridAverage', 0.0005);
figure;
axis equal;
hold on;
pcshow(rightWall);
pcshow(fixedRight, 'MarkerSize',400);
% return;
pcshow(leftWall);
pcshow(fixedLeft, 'MarkerSize',400);
return;


% ===================================================================================================

function [ptc, fu, fv, image] = GrabFrameData(data, frame)
    rgbW = 1280;
    rgbH = 720;
    depthW = 640;
    depthH = 480;
    rgb = data(frame,20);
    rgb = cell2mat(rgb.FrameColor);
    rgb = reshape(rgb,3,[]);
    r = reshape(rgb(1,:), rgbW, rgbH)'; 
    g = reshape(rgb(2,:), rgbW, rgbH)'; 
    b = reshape(rgb(3,:), rgbW, rgbH)';
    image = cat(3, r, g, b);
    %figure, imshow(image);
    x = (cell2mat(data(frame, 21).FrameX));
    y = (cell2mat(data(frame, 22).FrameY));
    z = (cell2mat(data(frame, 23).FrameZ));
    u = (cell2mat(data(frame, 24).FrameU));
    v = (cell2mat(data(frame, 25).FrameV));
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
        if (fuv(idx,1) > 0.1) & (fuv(idx,1) < 0.9) & (fuv(idx,2) > 0.1) & (fuv(idx,2) < 0.9)

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







