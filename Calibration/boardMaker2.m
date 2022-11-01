close all; clear all; clc;
squareSize = 200;
cb1 = checkerboard(squareSize,5,8);
cb1 = cb1(:,1:15*squareSize) ;
cb3 = (cb1 > 0.5);
r = uint8(255*ones(size(cb3)));
g = uint8(255*cb3);
b = uint8(255*cb3);
red = cat(3, r, g, b);
rx = randi([0,15], 20, 1);
ry = randi([0,3], 20, 1);
rr = randi([30,80], 20, 1);
figure, imshow(red);
hold on;
rad = 30;
halfSquareSize = squareSize * 0.5;
for i=1:10
    for j=1:15
        if (rem(i,2) == 0)
            if (rem(j,2) == 0)
                DrawCircles(squareSize, 10*squareSize- (squareSize*i - halfSquareSize), squareSize * j - halfSquareSize, rad, i );
            else
                DrawCircles(squareSize, 10*squareSize- (squareSize*i - halfSquareSize), squareSize * j - halfSquareSize, rad, j );
            end
        else
            if (rem(j,2) > 0)
                DrawCircles(squareSize, 10*squareSize- (squareSize*i - halfSquareSize), squareSize * j - halfSquareSize, rad, i );
            else
                DrawCircles(squareSize, 10*squareSize- (squareSize*i - halfSquareSize), squareSize * j - halfSquareSize, rad, j );
            end
        end
    end
end

% figure, imshow(cb1);
% hold on;
% for i=1:20
%     plot(rx(i)*100+50, ry(i)*100+50, '.r','MarkerSize', rr(i));
% end
exportgraphics(gcf,'left.png','Resolution',300);
F = getframe(gca);
[imBase Map] = frame2im(F);
close all;
[colorCorners1, boardSize1] = detectCheckerboardPoints(imBase);
tiledlayout(3,1);
nexttile;
imshow(imBase); hold on;
scatter(colorCorners1(:,1), colorCorners1(:,2),'oy', 'filled');

theta = 10;
tm = [cosd(theta) -sind(theta) 0.00050; ...
    sind(theta) cosd(theta) 0.0002; ...
    0 0 1];
tform = projective2d(tm);
imWarped = imwarp(imBase, tform);

nexttile; imshow(imWarped);  hold on;
[colorCorners2, boardSize2] = detectCheckerboardPoints(imWarped);
scatter(colorCorners2(:,1), colorCorners2(:,2),'oy', 'filled');

projectiveTransformation = fitgeotrans(colorCorners2, colorCorners1, 'Projective');
imFixed = imwarp(imWarped, projectiveTransformation);
stats = regionprops(imFixed, 'BoundingBox' );
T=struct2table(stats);

bbx = st.BoundingBox;
lx = floor(bbx(1));
ly = floor(bbx(2));
ux = lx + floor(bbx(4));
uy = ly + floor(bbx(5));
imFixed = imFixed(ly:uy, lx:ux, :);

nexttile; imshow(imFixed); hold on;
colorCorners2 = [colorCorners2;0,0]
[x,y] = transformPointsForward(projectiveTransformation, colorCorners2(:,1), colorCorners2(:,2));
colorCorners3 = [x,y ];scatter(colorCorners3(:,1), colorCorners3(:,2),'oy', 'filled');




cirb = (imFixed(:,:,1) > 0) & (imFixed(:,:,2) == 0);
[centers, radii] = imfindcircles(cirb,[5 50],'Method','TwoStage','Sensitivity',0.75)
viscircles(centers, radii,'EdgeColor','b');
[colorCorners3, boardSize] = detectCheckerboardPoints(imFixed);
scatter(colorCorners3(:,1), colorCorners3(:,2),'oy', 'filled');






function [xunit,yunit] = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end

function DrawCircles(sSize, y, x, r, n)
nx = ceil(sqrt(n));
ny = ceil( n / nx);
sx = sSize / (nx + 1);
sy = sSize / (ny + 1);
counter = 1;
for i=1:nx
    for j=1:ny
        if (counter <= n)
            plot(x - sSize * 0.5 + i * sx, y - sSize * 0.5 + j * sy, '.k','MarkerSize', r);
        end
        counter = counter + 1;
    end
end
end