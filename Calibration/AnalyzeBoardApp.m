close all;  clear all; clc;
image = imread("left2.bmp");
[squares, prevailingColor] = AnalyzeBoard2(image,1)
return;
image = imread("right2.bmp");
disp('right2');
[squares, prevailingColor] = AnalyzeBoard(image,1)
