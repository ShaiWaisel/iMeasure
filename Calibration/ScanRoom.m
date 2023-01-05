close all; clear all;imtool close all;
VERBOSE = 1;
platformID= "BaseP001";
%configFileName = 'C:\Projects\iMeasure\FromErez\Database\!iMeasure.db';
configFileName = 'C:\Projects\iMeasure\FromErez\Database\!iMeasure_LRB.db';
%dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SAMCore7114d.db';
%dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SAMCore_frameframe_13_11_22_4.db';
dataFileName = 'C:\Projects\iMeasure\FromErez\Database\ScanLRB180_1.db';
plyFileName = 'combined.ply';
SENSORS = 3;
FRAME = 1;

xlimits = [-2 2]; % meters
ylimits = [-2 8];
zlimits = [-2 2];




% Read sensors IDs and transformation matrices from configuration DB
connPlatform = sqlite(configFileName); 
C:\Projects\LightForce\ScanTrimmer
sqlQuery = sprintf('SELECT * FROM Platforms WHERE PlatformID = "%s"', platformID);
platformData=table2array(fetch(connPlatform, sqlQuery));
close(connPlatform);

% Allocate memory
sensors = cell(1, SENSORS);
transMats = cell(1, SENSORS);
rawData = cell(1, SENSORS);
wall = cell(1, SENSORS);
transWall = cell(1, SENSORS);

% Get data from DB structure
for i=1:SENSORS
    sensors{i} = platformData(5+i*4);
    transMats{i} = sscanf(platformData(8+i*4), '%f;', [4, inf])';
end

% Open data DB and read raw data
conn = sqlite(dataFileName); 
for sensor=1:length(sensors)
    rawData{sensor} = sqlread(conn,sprintf('%s_Frames',sensors{sensor}));
    wall{sensor}= [];
    transWall{sensor}= [];
end
close(conn);


FRAMES = size(rawData{1},1);
tform = affine3d;
loop = 0;
model = [];
fixed = pointCloud([0,0,0]);
fullFixed = fixed;
scannedModel = fixed;
goodShotIndex = 0;
for frame=1:13%FRAMES
    loop = loop + 1;

    clouds = [];
    blurred = 0;
    for sensor=1:SENSORS
        [cloud, u, v, img] = GrabFrameData(rawData{sensor}, loop, 0.1);
        mat = transMats{sensor}; 
        tf = affinetform3d(mat');
        cloud = pctransform(cloud, tf);
        frameCloud{sensor} = cloud;
        clouds = [clouds, cloud];
        blurred = max(blurred, blurMetric(img));
        if (blurred > 0.4)
            fprintf('Blurred frame %d Sensor %d\n', frame, sensor);
        end
    end
    original = pccat(clouds);
    [planes, modified, planars, bounds{loop}, edges{loop}, corners{loop}] = SetPlanarAreas(original, 0.02, 10000);
    if (size(planes,1) > 2) && (blurred < 0.4)
        goodShotIndex = goodShotIndex + 1;
        fprintf('Input frame  %d, points %d, planes %d\n',loop, size(original.Location,1), size(planes,1));
        fullMoving = pcdownsample(original, 'gridAverage', 0.01);
        model{goodShotIndex} = original;
        moving = pcdownsample(planars, 'gridAverage', 0.01);
        if (size(planes,1) > 2)
            if (size(fixed.Location,1)>1)
                tform = pcregistericp(moving, fixed, 'Metric','PlaneToPlane','Tolerance', [0.001 0.05] , 'InitialTransform', tform);
                %tform = pcregisterndt(moving, fixed, 0.2, 'InitialTransform', tform);
            end
        else
            if (size(fullFixed.Location,1)>1)
                %tform = pcregistericp(fullMoving, fullFixed, 'Metric','pointToPlane','Extrapolate', true, 'InitialTransform', tform);
                tform = pcregisterndt(fullMoving, fullFixed, 0.3, 'InitialTransform', tform);
            end
            
        end
        tMats{goodShotIndex} = tform; 
        for sensor=1:SENSORS
            tFrames{goodShotIndex, sensor} = frameCloud{sensor};
        end
        fixed = moving;
        fullFixed = fullMoving;
    end
end

accumTform = rigid3d;
scannedModel = model{1};

for frame = 2:size(tMats, 2)
    tform = tMats{frame};
    points = model{frame};
    accumTform.T = tform.T * accumTform.T;
    fixed = pcdownsample(scannedModel, 'gridAverage', 0.01);
    points = pctransform(points,accumTform);
    moving = pcdownsample(points, 'gridAverage', 0.01);
    tf = pcregistericp(moving, fixed, 'Metric','PlaneToPlane','Tolerance', [0.001 0.05] );
    %accumTform.T = tf.T * accumTform.T;
    points = pctransform(points,tf);
    %scannedModel = UnifyExtended(scannedModel, points);
    scannedModel = pccat([scannedModel, points]);
    fprintf('Accumulating frame  %d, points %d\n',frame, size(scannedModel.Location,1));

end
% scannedModel = pcdownsample(scannedModel, 'nonuniformGridSample', nonuniformGridSample);
% scannedModel.Normal = pcnormals(scannedModel); 
%pcshow(scannedModel,'Parent', lidarPlayer.Axes);


%show the result
player = pcplayer([-1.5 1.5], [-1.5 1.5], [-1 1]);
view(player.Axes, 10 , 80);

% Customize player axes labels
xlabel(player.Axes, 'X (m)')
ylabel(player.Axes, 'Y (m)')
zlabel(player.Axes, 'Z (m)')
hold(player.Axes, 'on');

view(player, scannedModel);
pcwrite(scannedModel,plyFileName);
return;


% ===================================================================================================

function [ptc, fu, fv, image] = GrabFrameData(data, frame, filterMargins)
    rgbW = 1280;
    rgbH = 720;
% build image raw data from DB 
    rgb = data(frame,19);
    rgb = cell2mat(rgb.FrameColor);
    rgb = reshape(rgb,3,[]);
    r = reshape(rgb(1,:), rgbW, rgbH)'; 
    g = reshape(rgb(2,:), rgbW, rgbH)'; 
    b = reshape(rgb(3,:), rgbW, rgbH)';
    image = cat(3, r, g, b);
    %figure, imshow(image);
%read XYZUV data from DB
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


function [planes, modified, planars, boundsList, edges, corners] = SetPlanarAreas(ptc, maxDistance, planeInPixels)
    referenceVector = [0,1,0];
    maxAngularDistance = 80;
    minColor = 100;
    inliers = ((ptc.Color(:,1) > minColor) & ...
                (ptc.Color(:,2) > minColor) & ...
                (ptc.Color(:,3) > minColor)) ;
    ptc = select(ptc, inliers);

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



