function [scannedModel] = StitchFrames(rawData, sensors, transMats, frameSpan, params)
    
    SENSORS = length(sensors);
%             params.PlaneError = app.editPlaneError.Value;
%             params.ICPgrid = app.editICPgrid.Value;
%             params.PointsTolerance = app.editPointsTolerance.Value;
%             params.Axes = {app.Sensor1Axis, app.Sensor2Axis, app.Sensor3Axis, app.Sensor4Axis};
%             params.Gauge = app.Gauge;
%             params.Margins
%             params.PlanePoints

    
    tform = affine3d;
    loop = 0;
    model = [];
    fixed.Location = [];
    fullFixed = fixed;
    goodShotIndex = 0;
    for frame=frameSpan(1):frameSpan(2)
        params.FrameGauge.Value = frame;
        drawnow;
        loop = loop + 1;
    
        clouds = [];
        blurred = 0;
        for sensor=1:SENSORS
            [cloud, u, v, img] = GrabFrameData(rawData{sensor}, loop, params.Margins);
            imshow(img, 'Parent', params.Axes{sensor});
            drawnow;
            frameCloud{sensor} = cloud;
            mat = transMats{sensor}; 
            tf = affinetform3d(mat');
            cloud = pctransform(cloud, tf);
            clouds = [clouds, cloud];
            blurred = max(blurred, blurMetric(img));
            if (blurred > 0.4)
                fprintf('Blurred frame %d Sensor %d\n', frame, sensor);
            end
        end
        black = img * 0;
            imshow(black, 'Parent', params.Axes{4});
        
        original = pccat(clouds);
    %    [planes, modified, planars, bounds{loop}, edges{loop}, corners{loop}] = SetPlanarAreas(original, 0.01, 20000);
        [planes, modified, planars, bounds{loop}, edges{loop}, corners{loop}] = SetPlanarAreas(original, params.PlaneError, params.PlanePoints, 0);
        if (size(planes,1) > 2) && (blurred < 0.4)
            goodShotIndex = goodShotIndex + 1;
            fprintf('Input frame  %d, points %d, planes %d\n',loop, size(original.Location,1), size(planes,1));
            fullMoving = pcdownsample(modified, 'gridAverage', params.ICPgrid);
            model{goodShotIndex} = modified;
            moving = pcdownsample(planars, 'gridAverage', params.ICPgrid);
            if (size(planes,1) > 2)
                if (size(fixed.Location,1)>1)
                    tform = pcregistericp(moving, fixed, 'Metric','PlaneToPlane','Tolerance', [params.PointsTolerance params.ICPgrid] , 'InitialTransform', tform);
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
    params.DisplayGauge.Value = 1;
    params.DisplayGauge.Visible = 1;
    params.DisplayGauge.Limits = [1, size(tMats,2)];
    
    for frame = 2:size(tMats, 2)
        tform = tMats{frame};
        tframe = tFrames{frame,:};
        scannedFrame = [];
        newFrame = [];
        for sensor=1:SENSORS
            newPoints = SelectExtended(tFrames{frame, sensor}, rotm2eul(tMats{frame}.T(1:3,1:3)));
            if (~isempty(newPoints.Location))
            %newPoints = pctransform(tFrames{frame, sensor}, tform);
                mat = transMats{sensor}; 
                tf = affinetform3d(mat');
                framePoints = pctransform(tFrames{frame, sensor}, tf);
                newPoints = pctransform(newPoints, tf);
                scannedFrame = pccat([scannedFrame, framePoints]);
                newFrame = pccat([newFrame, newPoints]);
            end
        end
        accumTform.T = tform.T * accumTform.T;
    
        scannedFrame = pctransform(scannedFrame, accumTform);
        points = pctransform(scannedFrame,accumTform);
        newFrame = pctransform(newFrame, accumTform);
    %     scannedModel = pccat([scannedModel, points]);
    %     accumTform.T = tform.T * accumTform.T;
    
        fixed = pcdownsample(scannedModel, 'gridAverage', params.ICPgrid);
        %points = pctransform(points,accumTform);
        moving = pcdownsample(scannedFrame, 'gridAverage', params.ICPgrid);
        tf = pcregistericp(moving, fixed, 'Metric','PlaneToPlane','Tolerance', [params.PointsTolerance params.ICPgrid] );
        %accumTform.T = tf.T * accumTform.T;
        newFrame = pctransform(newFrame, tf);
        scannedModel = pccat([scannedModel, newFrame]);
        params.DisplayGauge.Value = frame;
        drawnow;
    
    %     points = pctransform(points,tf);
    %         scannedModel = pccat([scannedModel, points]);
    
    
    
    
    %    points = pctransform(points,accumTform);
    %     scannedModel = UnifyExtended(scannedModel, points);
    %    scannedModel = pccat([scannedModel, points]);
        fprintf('Accumulating frame  %d, points %d\n',frame, size(scannedModel.Location,1));
    
    end
    params.DisplayGauge.Visible = 0;
    % scannedModel = pcdownsample(scannedModel, 'nonuniformGridSample', nonuniformGridSample);
    % scannedModel.Normal = pcnormals(scannedModel); 
    %pcshow(scannedModel,'Parent', lidarPlayer.Axes);
end


% ===================================================================================================

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


function [planes, modified, planars, boundsList, edges, corners] = SetPlanarAreas(ptc, maxDistance, planeInPixels, snapNormal)
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
            if (snapNormal) % project point to the closest point on the plane
                newLocs = planar.Location - plane.Normal .* distFromPlane;   
            else
                % project point along ray from sensor alongY
                points = planar.Location;
                pointsOnPlane = repmat([-plane.Parameters(4) /plane.Parameters(1) 0 0],size(points,1),1);
                rays = repmat([0 1 0],size(points,1),1);
                planeNormals = repmat([plane.Parameters(1) plane.Parameters(2) plane.Parameters(3)],size(points,1),1);
                % vector form of line_plane_intersection.m
                d1 = -dot(planeNormals,pointsOnPlane,2);
                t = -(d1 + dot(planeNormals,points,2)) ./ dot(planeNormals, rays,2);
                newLocs=points+rays.*t;
            end
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
%     for i=1:size(intersects,1)
%         P0 = intersects(i,1:3)-intersects(i,4:6) * span;
%         P1 = intersects(i,1:3)+intersects(i,4:6) * span;
%         edges = [edges; P0; P1;[NaN, NaN, NaN]];
%         for j=i+1:size(intersects,1)
%             i
%             j
%             Q0 = intersects(j,1:3)-intersects(j,4:6) * span;
%             Q1 = intersects(j,1:3)+intersects(j,4:6) * span;
%             [D,Xcp,Ycp,Zcp,Xcq,Ycq,Zcq,Dmin,imin,jmin]= ll_dist3d(P0, P1, Q0, Q1);
%             if (D < 0.01) && (Xcp > -span) && (Xcp < span) && (Ycp > -span) && (Ycp < span) && (Zcp > -span) && (Zcp < span)
%                 corners = [corners; Xcp, Ycp, Zcp];
%             end
%         end
%     end
end

function [selected] = SelectExtended(cloud, yawPitchRoll)
if (~isempty(cloud.Location))
    margins = [cloud.XLimits(1),cloud.XLimits(2),cloud.ZLimits(1),cloud.ZLimits(2) ];
    fovX = deg2rad(70);
    spanAng = asin(yawPitchRoll(1));
            
        if (yawPitchRoll(1) > 0)
            margins(2) = cloud.XLimits(2) - 2*spanAng/fovX* (cloud.XLimits(2) - cloud.XLimits(1));
        else
            margins(1) = cloud.XLimits(1) - 2*spanAng/fovX* (cloud.XLimits(2) - cloud.XLimits(1)); 
        end
        outliers = ((cloud.Location(:,1) < margins(1)) | ...
                    (cloud.Location(:,1) > margins(2)) | ...
                    (cloud.Location(:,3) < margins(3)) | ...
                    (cloud.Location(:,3) > margins(4))) ;
        selected = select(cloud, outliers);
else
    selected = cloud;
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



