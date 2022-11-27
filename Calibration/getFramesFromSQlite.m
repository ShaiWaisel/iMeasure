%close all; clear all;
%conn = sqlite('C:\Projects\iMeasure\FromErez\Database\SAMCor-30Calib.db'); data = sqlread(conn,'f1230591_Frames');
%data = sqlread(conn,'f0171826_Frames');

conn = sqlite('C:\Projects\iMeasure\FromErez\Database\RedCircles\SAMCore_red_circles_first_scan.db'); data = sqlread(conn,'f0171826_Frames');
close(conn);
FRAMES = min(size(data,1),15);

xlimits = [-2 2]; % meters
ylimits = [-2 8];
zlimits = [-2 2];

lidarPlayer = pcplayer(xlimits, ylimits, zlimits,'VerticalAxis', 'Z');
view(lidarPlayer.Axes, 10 , 80);

% Customize player axes labels
xlabel(lidarPlayer.Axes, 'X (m)')
ylabel(lidarPlayer.Axes, 'Y (m)')
zlabel(lidarPlayer.Axes, 'Z (m)')



gridSize = 0.01;
nonuniformGridSample = 6;
h2=[];
fixed = pointCloud([0,0,0]);
fullFixed = fixed;
scannedModel = fixed;
accumTform = affine3d;
Nframes = 41;
transMats = affine3d;
model = [];
tform = affine3d;
loop = 0;
for frame=1:1:FRAMES
    loop = loop + 1;
    [original, fu, fv, image] = GrabFrameData(data, frame);
    [planes, modified, planars, bounds{loop}, edges{loop}, corners{loop}] = SetPlanarAreas(original, 0.02, 20000);
    fullMoving = pcdownsample(original, 'gridAverage', 0.005);
    model{loop} = fullMoving;
    moving = pcdownsample(planars, 'gridAverage', 0.01);
    fprintf('Loop %d, planes %d\n',loop, size(planes,1));
    if (size(planes,1) > 2)
        if (size(fixed.Location,1)>1)
            tform = pcregistericp(moving, fixed, 'Metric','pointToPlane','Extrapolate', true, 'InitialTransform', tform);
            %tform = pcregisterndt(moving, fixed, 0.3, 'InitialTransform', tform);
        end
    else
        if (size(fullFixed.Location,1)>1)
            %tform = pcregistericp(fullMoving, fullFixed, 'Metric','pointToPlane','Extrapolate', true, 'InitialTransform', tform);
            tform = pcregisterndt(fullMoving, fullFixed, 0.3, 'InitialTransform', tform);
        end
        
    end
    tMats{loop} = tform; 
    fixed = moving;
    fullFixed = fullMoving;
end

hold(lidarPlayer.Axes, 'on');
accumTform = rigid3d;
scannedModel = model{1};
for frame = 2:size(tMats, 2)
    tform = tMats{frame};
    points = model{frame};
    accumTform.T = tform.T * accumTform.T;
    points = pctransform(points,accumTform);
    scannedModel = UnifyExtended(scannedModel, points);
    fprintf('Accumulating frame  %d, points %d\n',frame, size(scannedModel.Location,1));
    %pcshow(points,'Parent', lidarPlayer.Axes);
   
%         if (size(h2,1) > 0)
%         for i=1:size(h2,1)
%             delete (h2(i));
%         end
%         h2 = [];
%     end;
%     for i=1:size(bounds,2)
%         bounds = bounds{i};
%         h2 = [h2; plot3(bounds(:,1), bounds(:,2), bounds(:,3),'Parent',lidarPlayer.Axes,'Color',[1.0, 1.0, 0.0])];
%     end


end
% scannedModel = pcdownsample(scannedModel, 'nonuniformGridSample', nonuniformGridSample);
% scannedModel.Normal = pcnormals(scannedModel); 
pcshow(scannedModel,'Parent', lidarPlayer.Axes);
return;
hold(lidarPlayer.Axes, 'on');
[modified, bounds, edges, corners] = SetPlanarAreas(scannedModel, 0.05);
scatter3(lidarPlayer.Axes,corners(:,1), corners(:,2), corners(:,3)); plot3(lidarPlayer.Axes,edges(:,1), edges(:,2), edges(:,3))
pcshow(modified,'Parent', lidarPlayer.Axes);

transMats


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
    indices = Features(fu, fv, image, [778, 295]);
end

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

function [edges, corners] = transformEdgesCorners(e, c, tform)
e = [e, ones(sizeof(e,1))];
edges = e * tform.T;
c = [c, ones(sizeof(c,1))];
corners = c * tform.T;
end



function [bounds] = PlaneBounds(pts, pe)
    idx = randi([1, size(pts,1)],1000,1);
    b = boundary(pts(idx,1), pts(idx,3));
    ib = idx(b);
    bpy = -(pe(1)*pts(ib,1) + pe(3)*pts(ib,3) + pe(4))/pe(2);
    bounds = [pts(ib,1), bpy, pts(ib,3)];
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

function [featuresIndices] = Features(cloudU, cloudV, image, featuersPos)

    fuv = [featuersPos(:,1) / size(image,2), featuersPos(:,2) / size(image,1)];
    featuresIndices = [];
    for idx=1:size(fuv,1)
        dists = (fuv(idx, 1) - cloudU).*(fuv(idx,1) -cloudU) + (fuv(idx, 2) - cloudV).*(fuv(idx,2) -cloudV);
        
        [d,idx] = min(dists);
        featuresIndices = [featuresIndices, idx];
    end
end






