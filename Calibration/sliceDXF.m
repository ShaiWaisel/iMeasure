function [poly] = sliceDXF(ptc, Z, tolerance, dxfFileName)
stepx = (ptc.XLimits(2) - ptc.XLimits(1)) / 100;
stepy = (ptc.YLimits(2) - ptc.YLimits(1)) / 100;
poly = [];
neighbors = 5;
for x=ptc.XLimits(1):stepx:ptc.XLimits(2)
    testP = [x,ptc.YLimits(1),Z];
    [TestPointsIndices, dists] = findNearestNeighbors(ptc, testP,neighbors); %averaging K closest neighbours
    m = mean(ptc.Location(TestPointsIndices, :),1);
    if (norm(m-testP) < tolerance)
        poly = [poly; m(1:2),Z];
    end
end
for y=ptc.YLimits(1):stepy:ptc.YLimits(2)
    testP = [ptc.XLimits(2), y, Z];
    [TestPointsIndices, dists] = findNearestNeighbors(ptc, testP,neighbors); %averaging K closest neighbours
    m = mean(ptc.Location(TestPointsIndices, :),1);
    if (norm(m-testP) < tolerance)
        poly = [poly; m(1:2),Z];
    end
end
for x=ptc.XLimits(2):-stepx:ptc.XLimits(1)
    testP = [x,ptc.YLimits(2),Z];
    [TestPointsIndices, dists] = findNearestNeighbors(ptc, testP,neighbors); %averaging K closest neighbours
    m = mean(ptc.Location(TestPointsIndices, :),1);
    if (norm(m-testP) < tolerance )
        poly = [poly; m(1:2),Z]; 
    end
end
for y=ptc.YLimits(2):-stepy:ptc.YLimits(1)
    testP = [ptc.XLimits(1), y, Z];
    [TestPointsIndices, dists] = findNearestNeighbors(ptc, testP,neighbors); %averaging K closest neighbours
    m = mean(ptc.Location(TestPointsIndices, :),1);
    if (norm(m-testP) < tolerance)
        poly = [poly; m(1:2),Z];
    end
end

if (size(poly,1) > 0)
    FID = dxf_open(dxfFileName);
    dxf_polyline(FID,poly(:,1), poly(:,2), repmat([Z],size(poly,1),1));
    dxf_close(FID);
end
end