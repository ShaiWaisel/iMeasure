function [ptcAligned, limits] = AlignCloudToAxes(ptc, VERBOSE)
    maxAngularDistance = 45;
    maxDistance = 0.03;
    refAxis = [1 0 0];
    [plane,inlierIndices,outlierIndices] = pcfitplane(ptc, maxDistance,refAxis,maxAngularDistance);
    p1 = select(ptc, inlierIndices);
    if (VERBOSE)
        figure; pcshow(p1);axis equal;
    end
    ax = cross(plane.Normal, refAxis);
    ang = acos(dot(refAxis, plane.Normal));
    axisa = [ax, ang];
    tf = axang2tform(axisa);
    tm = affinetform3d(tf);
    ptc = pctransform(ptc,tm);
    p1 = pctransform(p1,tm);
    refAxis = [0 1 0];
    [plane,inlierIndices,outlierIndices] = pcfitplane(ptc, maxDistance,refAxis,maxAngularDistance);
    p2 = select(ptc, inlierIndices);
    if (VERBOSE)
        figure; pcshow(p2);axis equal;
    end
    ax = cross(plane.Normal, refAxis);
    ang = acos(dot(refAxis, plane.Normal));
    axisa = [ax, ang];
    tf = axang2tform(axisa);
    tm = affinetform3d(tf);
    ptcAligned = pctransform(ptc,tm);
    refAxis = [1 0 0];
    planes2go = 1;
    planeInPixels = 50000;
    planars = [];
    maxAngularDistance = 5;
    ptc = ptcAligned;
    while (planes2go > 0)
        planes2go = 0;
        [plane,inlierIndices,outlierIndices] = pcfitplane(ptc, maxDistance,refAxis,maxAngularDistance);
        if (~isempty(inlierIndices))
            if (isempty(planars))
                planars = ptc.select(inlierIndices);
            else
                planars = pccat([planars, ptc.select(inlierIndices)]);
            end
            ptc = ptc.select(outlierIndices);
            if (size(inlierIndices, 1) > planeInPixels)
                planes2go = 1;
            end
        end
    end
    refAxis = [0 1 0];
    planes2go = 1;
    ptc = ptcAligned;
    while (planes2go > 0)
        planes2go = 0;
        [plane,inlierIndices,outlierIndices] = pcfitplane(ptc, maxDistance,refAxis,maxAngularDistance);
        if (~isempty(inlierIndices))
            if (isempty(planars))
                planars = ptc.select(inlierIndices);
            else
                planars = pccat([planars, ptc.select(inlierIndices)]);
            end
            ptc = ptc.select(outlierIndices);
            if (size(inlierIndices, 1) > planeInPixels)
                planes2go = 1;
            end
        end
    end

    limits = [planars.XLimits; planars.YLimits; planars.ZLimits];
    tf = eye(4);
    tf(1:3,4) = -limits(:,1)'; 
    tm = affinetform3d(tf);
    ptcAligned = pctransform(ptcAligned,tm);
    if (VERBOSE)
        pp=pcplayer([-5 5],[-5 5],[-2 2]);
        view(pp,p12);
    end
end

