% Analyze board with prevailing colors squares and black circles
function [squares, prevailingColor] = AnalyzeBoard(image, verbose)
[xnd,map] = rgb2ind(image,3,'nodither');
img3Colors = ind2rgb(xnd, map);
R = sum(sum(img3Colors(:,:,1)));
G = sum(sum(img3Colors(:,:,2)));
B = sum(sum(img3Colors(:,:,3)));
M = max([R,G,B]);
N = min([R,G,B]);
imPrevail = img3Colors;
switch (M)
    case R
        prevailingColor = 1;
        imPrevail = (img3Colors(:,:,1)>img3Colors(:,:,2)*2) & (img3Colors(:,:,1)>img3Colors(:,:,3)*2) & (img3Colors(:,:,1) > 0.4);
    case G
        prevailingColor = 2;
        imPrevail = (img3Colors(:,:,2)>img3Colors(:,:,1)*2) & (img3Colors(:,:,2)>img3Colors(:,:,3)*2) & (img3Colors(:,:,2) > 0.4);
    case B
        prevailingColor = 3;
        imPrevail = (img3Colors(:,:,3)>img3Colors(:,:,1)*2) & (img3Colors(:,:,3)>img3Colors(:,:,2)*2) & (img3Colors(:,:,3) > 0.4);
    otherwise
        prevailingColor = 4;
end

prevailingAOI = imerode(imPrevail,strel('disk',7));
stats = regionprops(prevailingAOI,"BoundingBox");
T=struct2table(stats);
minX = round(max(1,min(T.BoundingBox(:,1))));
maxX = round(min(size(image,2),max(T.BoundingBox(:,1)+T.BoundingBox(:,3))));
img3Colors = img3Colors(:,minX:maxX,:);

imgGray = uint8(rgb2gray(img3Colors)*255);
[corners, boardSize] = detectCheckerboardPoints(imgGray, 'HighDistortion',true,'PartialDetections',true, 'MinCornerMetric',0.18);
x = reshape(corners(:,1),boardSize(1)-1,[]);
y = reshape(corners(:,2),boardSize(1)-1,[]);
if (verbose)
    figure, imshow(img3Colors);hold on;
    figure, imshow(image);hold on;
end
NCols = size(x,2);
NRows = size(x,1);
square = struct('color',0,'value',0, 'quad', zeros(4,2),'row',0, 'col', 0);
[squares(1:NRows, 1:NCols)] =deal(square);
%% Set color and value per square
for row=1:NRows-1
    for col=1:NCols-1
    quad = [x(row,col)+1, y(row,col)+1; x(row+1, col)+1, y(row+1, col); x(row+1, col+1), y(row+1, col+1); x(row, col+1), y(row, col+1)+1];
    scatter(quad(:,1)+minX, quad(:,2),50,'oy', 'filled');
    if (verbose)
        [color, value] = GetSquareValues(img3Colors, quad,0);
    end
    fprintf('col=%d row=%d color=%d value=%d\n', col, row, color, value);
    squares(row,col).color =color;
    squares(row,col).value =value;
    squares(row,col).quad = quad+ minX;
    end
end
%% Set Row value as per the whites
for row=1:NRows-1
    whites=[];
    for col=1:NCols-1
        if(squares(row,col).color == 4)
            whites = [whites; squares(row,col).value];
        end
    end
    if (size(whites,1)>0)
        rowNumber = round(mean(whites));
        for col=1:NCols-1
            squares(row,col).row = rowNumber;
        end
    end
end

%% Set Col value as per the colored
for col=1:NCols-1
    colored=[];
    for row=1:NRows-1
        if(squares(row,col).color ~= 4)
            colored = [colored; squares(row,col).value];
        end
    end
    if (size(colored,1)>0)
        colNumber = round(mean(colored));
        for row=1:NRows-1
            squares(row,col).col = colNumber;
        end
    end
end

end


function [color, value] = GetSquareValues(image, square, verbose)
x = square(:,1)'; y = square(:,2)';
x = [x, x(1)];
y = [y, y(1)];
[rows, columns, numberOfColorChannels] = size(image);
mask = poly2mask(x, y, rows, columns);
mask3 = cat(3, ~mask, ~mask, ~mask);
image(mask3) = 0;
x1 = floor(min(x));
x2 = floor(max(x));
y1 = floor(min(y));
y2 = floor(max(y));
cropped = image(y1:y2,x1:x2,:);
% mask = poly2mask([x - min(x)+1], [y - min(y)+1], floor(max(y)-min(y)+1), floor(max(x)-min(x)+1));
% cropped(~mask,:) = 0;

if (verbose)
    figure, imshow(cropped);
end
w = 128;
h = 128;
r2 = [0, 0; 0, size(cropped,2); size(cropped,1), size(cropped,2); size(cropped,1), 0];
r2 = [0, 0; 0, w; h, w; h, 0];
sq = [square(:,1) - min(square(:,1)), square(:,2) - min(square(:,2))];
tf = fitgeotrans(sq, r2, 'projective');
ci = imwarp(cropped,tf);
if (verbose)
    figure, imshow(ci);
end
R = sum(sum(ci(:,:,1)));
G = sum(sum(ci(:,:,2)));
B = sum(sum(ci(:,:,3)));
M = max([R,G,B]);
N = min([R,G,B]);
if ((R/M == 1) && (R>N*1.5))
    color = 1;
    ci = ci(:,:,1);
else
    if ((G/M == 1) && (G>N*1.5))
        color = 2;
        ci = ci(:,:,2);
    else
        if ((B/M == 1) && (B>N*1.5))
            color = 3;
            ci = ci(:,:,3);
        else
            color = 4;
            ci = rgb2gray(ci);
        end
    end
end
ci = 1 - ci;
ci(ci<0.6) = 0;
ci(ci>0.6) = 1.0;
bi = im2bw(ci);
dotAreas = regionprops('table', bi,'Area');
da = (dotAreas.Area>100) & (dotAreas.Area<500);
value = sum(da);
end

