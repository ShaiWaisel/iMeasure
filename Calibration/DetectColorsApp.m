function DetectColors()
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.

image = imread('1left.bmp');
DetectColorsInImage(image, 1)
image = imread('1right.bmp');
DetectColorsInImage(image, 2)
image = imread('1back.bmp');
DetectColorsInImage(image, 3)
end

function DetectColorsInImage(image, row)
    columns = 7;
    rows = 3;
    R=image(:,:,1);
    G=image(:,:,2);
    B=image(:,:,3);
    g=rgb2gray(image);
    r1=imsubtract(R,g);
    g1=imsubtract(G,g);
    b1=imsubtract(B,g);

    y1 = (r1>5)  & (g1 > 5);
    rm = max(max(r1));
    imsubtract(r1, uint8(y1).*rm);
    rg = max(max(g1));
    imsubtract(g1, uint8(y1).*rg);
    r1 = (r1>30);
    g1 = (g1>20);
    b1 = (b1>20);
    r1 = bwareaopen(r1, 100);
    g1 = bwareaopen(g1, 100);
    b1 = bwareaopen(b1, 100);
    subplot(rows,columns,(row-1)*columns+1);
    imshow(image);
    subplot(rows,columns,(row-1)*columns+2);
    imshow(r1,[]);
    subplot(rows,columns,(row-1)*columns+3);
    imshow(g1,[]);
    subplot(rows,columns,(row-1)*columns+4);
    imshow(b1,[]);
    subplot(rows,columns,(row-1)*columns+5);
    imshow(y1,[]);
    col=5;
    s=sprintf('R=%d   G=%d    B=%d    Y=%d\n', sum(sum(r1)), sum(sum(g1)), sum(sum(b1)), sum(sum(y1)));
    if (sum(sum(r1)) > 10000)
        col=col+1;
        imageR = Crop(image, r1);
        subplot(rows,columns,(row-1)*columns+col);
        imshow(imageR,[]);
    end

    if (sum(sum(g1)) > 10000)
        col=col+1;
        imageG = Crop(image, g1);
        subplot(rows,columns,(row-1)*columns+col);
        imshow(imageG,[]);
    end

    if (sum(sum(b1)) > 10000)
        col=col+1;
        imageB = Crop(image, b1);
        subplot(rows,columns,(row-1)*columns+col);
        imshow(imageB,[]);
    end

    if (sum(sum(y1)) > 10000)
        col=col+1;
        imageY = Crop(image, y1);
        subplot(rows,columns,(row-1)*columns+col);
        imshow(imageY,[]);
    end
disp(s);
    
end

function [croppedImage] = Crop(image, mask)
p  = regionprops (mask, 'BoundingBox');
a = struct2array(p);
a = reshape(a,4,[])';
MIN=min(a(:,1:2));
r=[a(:,1)+a(:,3), a(:,2)+a(:,4)];
MAX=max(r);
x1 = ceil(MIN(1));
x2 = floor(MAX(1));
y1 = ceil(MIN(2));
y2 = floor(MAX(2));
croppedImage = image(ceil(MIN(2)):floor(MAX(2)),ceil(MIN(1)):floor(MAX(1)),:);
%figure; imshow(r1); hold on;rectangle('Position',[MIN,MAX-MIN],'EdgeColor','y');

end