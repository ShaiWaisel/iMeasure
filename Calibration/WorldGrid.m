function [rGrid, gGrid, bGrid, yGrid] = WorldGrid(squareSize, boardSize, cameraPos, lookAt, verbose)
    grid = meshgrid(1:15,1:10,1:3);
    for i=1:10
        for j=1:15
            grid(i,j,1) = squareSize*j - squareSize / 2;
            grid(i,j,2) = 0;
            grid(i,j,3) = squareSize*i - squareSize / 2;
        end
    end
    
    yGrid = Transform(grid, pi, [boardSize, 0] , cameraPos, lookAt);
    rGrid = Transform(grid, 0, [0,boardSize], cameraPos, lookAt);
    gGrid = Transform(grid, -pi/2,[boardSize,boardSize], cameraPos, lookAt);
    bGrid = Transform(grid, pi/2,[0,0], cameraPos, lookAt);
    if (verbose)
        figure; axis equal; hold on;
        scatter3(rGrid(:,:,1), rGrid(:,:,2), rGrid(:,:,3),'r', 'filled');
        scatter3(gGrid(:,:,1), gGrid(:,:,2), gGrid(:,:,3),'g', 'filled');
        scatter3(bGrid(:,:,1), bGrid(:,:,2), bGrid(:,:,3),'b', 'filled');
        scatter3(yGrid(:,:,1), yGrid(:,:,2), yGrid(:,:,3),'y', 'filled');
    end
end

function [transformed] = Transform(grid, rotWall, transWall, cameraPos, lookAt)
% transform wall to it's place on calibration layout
axang = [0 0 1 rotWall];
mat = eye(4);
mat(1:3,1:3) = axang2rotm(axang);
mat(1:2,4) = transWall;

axang = [0, 0, 1, pi/2-atan2(lookAt(2) - cameraPos(2), lookAt(1) - cameraPos(1))];
mat2 = eye(4);
mat2(1:3,1:3) = axang2rotm(axang);

for (i=1:size(grid,1))
    for (j=1:size(grid,2))
        point = [grid(i,j,1), grid(i,j,2), grid(i,j,3),  1];
        point = mat * point';
        point(1:3) = point(1:3) - cameraPos';
        point =  mat2 * point ;
        transformed(i,j,:) = point(1:3);
    end
end
end