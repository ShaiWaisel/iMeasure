% Analyze board with black squares and prevailing colors circles
function [squares, prevailingColor] = AnalyzeBoard(image, verbose)
[xnd,map] = rgb2ind(image,3,'nodither');
img3Colors = ind2rgb(xnd, map);
R = sum(sum(img3Colors(:,:,1)));
G = sum(sum(img3Colors(:,:,2)));
B = sum(sum(img3Colors(:,:,3)));
M = max([R,G,B]);
N = min([R,G,B]);
[corners, boardSize] = detectCheckerboardPoints(img3Colors, 'HighDistortion',true,'PartialDetections',true, 'MinCornerMetric',0.18);
minX = round(max(1,min(corners(:,1)-3)));
maxX = round(min(size(image,2),max(corners(:,1)+3)));
img3Colors = img3Colors(:,minX:maxX,:);
M=R;

switch (M)
    case R
        prevailingColor = 1;
    case G
        prevailingColor = 2;
    case B
        prevailingColor = 3;
    otherwise
        prevailingColor = 4;
end

x = reshape(corners(:,1)-minX,boardSize(1)-1,[]);
y = reshape(corners(:,2),boardSize(1)-1,[]);
verbose = 1;
if (verbose)
    figure, imshow(img3Colors);hold on;
    scatter(x,y,'yo', 'filled');
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
        if (sum(any(isnan(quad))) == 0)
            scatter(quad(:,1)+minX, quad(:,2),50,'oy', 'filled');
            if (verbose)
                value = GetSquareValues(img3Colors, quad,prevailingColor, 0);
            end
            fprintf('col=%d row=%d color=%d value=%d\n', col, row, prevailingColor, value);
            squares(row,col).color =prevailingColor;
            squares(row,col).value =value;
            squares(row,col).quad = quad+ minX;
        end
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


function [value] = GetSquareValues(image, square, prevailingColor, verbose)
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
w = 96;
h = 96;
r2 = [0, 0; 0, size(cropped,2); size(cropped,1), size(cropped,2); size(cropped,1), 0];
r2 = [0, 0; 0, w; h, w; h, 0];
sq = [square(:,1) - min(square(:,1)), square(:,2) - min(square(:,2))];
tf = fitgeotrans(sq, r2, 'projective');
ci = imwarp(cropped,tf, OutputView = imref2d(size(image))); 
if (verbose)
    figure, imshow(ci);
end
switch (prevailingColor)
    case 1
        bi = (ci(:,:,1) > 0.5) & (ci(:,:,2)<0.2) & (ci(:,:,3) < 0.2);        
    case 2
        bi = (ci(:,:,1) < 0.2) & (ci(:,:,2)>0.5) & (ci(:,:,3) < 0.2);        
    case 3
        bi = (ci(:,:,1) < 0.2) & (ci(:,:,2)<0.2) & (ci(:,:,3) > 0.5);        
    otherwise
        bi = rgb2bw(ci);
end

% ci(ci<0.6) = 0;
% ci(ci>0.6) = 1.0;
% bi = im2bw(ci);
rad = round((x2-x1)/16);
se = strel('diamond',rad);
sd = strel('diamond', rad-1);
ei = imerode(bi,se);
di = imdilate(ei,sd);
dotAreas = regionprops('table', di,'Area');
da = (dotAreas.Area>7);
value = sum(da);
end

