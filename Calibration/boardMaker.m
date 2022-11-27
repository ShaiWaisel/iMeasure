close all; clear all; clc;
cb1 = checkerboard(100,2,4);
cb3 = [cb1, cb1(:,1:700)];
s3 = cb3(:, 401:800);
 s3(s3 < 0.5) = 1.0;
 %s3(s3 < 0.5) = 0.7;
% s3(s3<1.0) = 0.7;
%cb3(:,401:800) = s3;
cb3(:,801:1200) = s3 ;
s3(s3 < 0.8) = 0.0;
cb3(:,1201:1500) = s3(:,1:300) ;
rx = randi([0,15], 20, 1);
ry = randi([0,3], 20, 1);
rr = randi([30,80], 20, 1);
figure, imshow(cb3);
hold on;
for i=1:20
    plot(rx(i)*100+50, ry(i)*100+50, '.r','MarkerSize', rr(i));
end
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
cirb = (imWarped(:,:,1) > 0) & (imWarped(:,:,2) == 0);
[centers, radii] = imfindcircles(cirb,[5 50],'Method','TwoStage','Sensitivity',0.75)
viscircles(centers, radii,'EdgeColor','b');

projectiveTransformation = fitgeotrans(colorCorners2, colorCorners1, 'Projective');
imFixed = imwarp(imWarped, projectiveTransformation);
st = regionprops(imFixed, 'BoundingBox' );
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