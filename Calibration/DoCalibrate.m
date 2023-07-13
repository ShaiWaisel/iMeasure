function [paramCode, transMats, transWalls] = DoCalibrate(scanFileName, paramsFileName, outputFileName, VERBOSE, app)
if (~exist(paramsFileName, 'file'))
    paramCode = -1;
    return;
end
rawData = 0;

[paramCode, sensors, devices, locations, transMats, rawData] = ReadFromSQLscan(scanFileName);
if (paramCode<0)
    return;
end
fid=fopen(paramsFileName);
C = textscan(fid,'%f');
fclose(fid);


for sensor=1:length(sensors)
    wall{sensor}= [];
    transWalls{sensor}= [];
    Uval{sensor} = [];
    Vval{sensor} = [];
    moving{sensor} = [];
    fixed{sensor} = [];
    transMats{sensor} = [];
    if (app ~= 0)
        app.TabGroup.Children(sensor).Title =sensors{sensor};
    end
end

% set global calibration
sensorPosition = C{1}(1:3)';
lookAtPoint = C{1}(4:6)';
squareSize = C{1}(7);
marginsSize = C{1}(8);
boardSize = 1.5;
if (app ~= 0)
    app.editCameraX.Value = sensorPosition(1);
    app.editCameraY.Value = sensorPosition(2);
    app.editCameraZ.Value = sensorPosition(3);
    app.editLookAtX.Value = lookAtPoint(1);
    app.editLookAtY.Value = lookAtPoint(2);
    app.editLookAtZ.Value = lookAtPoint(3);
    app.editSquareSize.Value = squareSize;
    app.editFrameMargins.Value= marginsSize * 100.0;
end


%create theoretical grid (centers of squares)
[worldGrid.rGrid, worldGrid.gGrid, worldGrid.bGrid, worldGrid.yGrid] = ...
    WorldGrid(squareSize, boardSize, sensorPosition, lookAtPoint, VERBOSE);

FRAMES = min(8,size(rawData{1}.ID,1));
tform = affine3d;
first = 0;
ax = 0;
for frame=1:FRAMES
    blurred = 0;
    for sensor=1:length(sensors)
        [original{sensor}, u{sensor}, v{sensor}, img{sensor}] = GrabFrameData(rawData{sensor}, frame, marginsSize);
        blurred = max(blurred, blurMetric(img{sensor}));
    end
    if (blurred > 0.38)
         fprintf('Blurred frame %d Sensor %d\n', frame, sensor);
    else
        for sensor=1:length(sensors)
            if (~first)
                wall{sensor} = original{sensor};
                Uval{sensor} = u{sensor};
                Vval{sensor} = v{sensor};
            else
                wall{sensor} = pccat([wall{sensor}, original{sensor}]);
                Uval{sensor} = [Uval{sensor}; u{sensor}];
                Vval{sensor} = [Vval{sensor}; v{sensor}];
            end
        end
        first = frame;
    end
end
if (VERBOSE)
    figure; axis equal; hold on; 
    hAxes = gca;
end
for sensor=1:length(sensors)
    if (app ~= 0)
        app.TabGroup.SelectedTab = app.TabGroup.Children(sensor);
        ax = findobj(app.TabGroup.SelectedTab.Children,'Type','axes');
        table = findobj(app.TabGroup.SelectedTab.Children,'Type','uitable');
        table.Data = eye(4);
        drawnow;
    end

    [~, ~, ~, image{sensor}] = GrabFrameData(rawData{sensor},first, marginsSize);
    if (app ~= 0)
        imshow(image{sensor},'Parent', ax);
        hold(ax,'on');
    end
    r={};
    e={};
    MIN_PIX_SPOT = 50;
    for colorCode=1:4
        if (size(transMats{sensor},1) >= 0)
            [code, ROI, mask] = DetectColors(image{sensor},colorCode, MIN_PIX_SPOT, VERBOSE);
            if (code == 1)
                % found colorCode area
                [squares, code] =  AnalyzeBoard3(image{sensor}, ROI, mask, colorCode, ax, VERBOSE);
                if (code == 1)
                    for i=1:size(squares,1)
                        for j=1:size(squares,2)
                            if (squares(i,j).value > 0)
                                center2D = mean(squares(i,j).quad);
                                featuresIdxs = Features(Uval{sensor}, Vval{sensor}, image{sensor}, center2D);
                                if (size(featuresIdxs,1) > 0)
                                    center3D = wall{sensor}.Location(featuresIdxs,:);
                                    moving{sensor} = [moving{sensor}; center3D];
                                    switch colorCode
                                        case 1
                                            fixed{sensor} = [fixed{sensor}; reshape(worldGrid.rGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                        case 2 
                                            fixed{sensor} = [fixed{sensor}; reshape(worldGrid.gGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                        case 3 
                                            fixed{sensor} = [fixed{sensor}; reshape(worldGrid.bGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                        case 4 
                                            fixed{sensor} = [fixed{sensor}; reshape(worldGrid.yGrid(squares(i,j).row, squares(i,j).col, :),[], 3)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    m=moving{sensor};
    f=fixed{sensor};
    if (size(m,1)>2)
        [mat, r{sensor}, e{sensor}] = absor(m', f');
        %f1=[f, ones(size(f,1),1)];
        m1=[m, ones(size(m,1),1)];
        transMats{sensor} = mat.M';
        z=m1*mat.M';
        if (VERBOSE)
            scatter3(ax, f(:,1), f(:,2), f(:,3),'ob'); scatter3(ax, m(:,1), m(:,2), m(:,3),'+r');
            for i=1:size(f,1)
                text(hAxes, f(i,1),f(i,2),f(i,3)+50,num2str(i));
                    text(sz, m(i,1),m(i,2),m(i,3)+50,num2str(i));
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
if (app ~= 0) && (length(e)>0)
    app.lblMat1err.Text = sprintf('Error (m):   %f',e{1}.errmax);
    app.lblMat1det.Text = sprintf('Determinant: %f',det(transMats{1}));
    app.lblMat2err.Text = sprintf('Error (m):   %f',e{2}.errmax);
    app.lblMat2det.Text = sprintf('Determinant: %f',det(transMats{2}));
    app.lblMat3err.Text = sprintf('Error (m):   %f',e{3}.errmax);
    app.lblMat3det.Text = sprintf('Determinant: %f',det(transMats{3}));
    drawnow;
end
failed = 0;
for sensor=1:length(sensors)
    if (size(transMats{sensor},1) == 0)
        failed = 1;
        paramCode = -1;
    end
end
fid = fopen(outputFileName,'wt');
if (failed)
fprintf(fid, 'Could not detect checkerboard calibration,0,Failed\n');
else    
fprintf(fid, 'Completed,100,OK\n');
end;
if (~failed)
    for sensor=1:length(sensors)
        mat = transMats{sensor}; 
        numbers = reshape(mat,1,[]); 
        line = sprintf("%s,%s,%s,16,%s\n",sensors{sensor},devices{sensor},locations{sensor},outMat(numbers));
        fprintf(fid, line);
        tform = affinetform3d(mat');
        transWalls{sensor} = pctransform(wall{sensor}, tform);
    end
end
fclose(fid);
end

% ===================================================================================================

function str = outMat(vec)
str = sprintf(';%1.5f',vec);
str = str(2:end);
end

function [ptc, fu, fv, image] = GrabFrameData(data, frame, filterMargins)
    rgbW = 1280;
    rgbH = 720;
% build image raw data from DB 
    rgb = data(frame,20);
    rgb = cell2mat(rgb.FrameColor);
    rgb = reshape(rgb,3,[]);
    r = reshape(rgb(1,:), rgbW, rgbH)'; 
    g = reshape(rgb(2,:), rgbW, rgbH)'; 
    b = reshape(rgb(3,:), rgbW, rgbH)';
    image = cat(3, r, g, b);
    %figure, imshow(image);
%read XYZUV data from DB
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
%filter out-of-image UV indices    
    idx = (fu<filterMargins) | (fu>(1.0-filterMargins)) | (fv<filterMargins) | (fv>(1.0-filterMargins));
    fx(idx)=[];
    fy(idx)=[];
    fz(idx) = [];
    fu(idx) = [];
    fv(idx) = [];
    iu = int32(floor(size(image,2)*fu)+1);
    iv = min(size(image,1),int32(floor(size(image,1)*fv)+1));
    linIdx = sub2ind(size(image), iv, iu);
%create points cloud in platform coordinates Z <= -Y, Y <= Z    
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







