% Analyze board with black squares and prevailing colors circles
function [squares, code] = AnalyzeBoard3(bigImage, ROI, imageBW, prevailingColor, ax, verbose)
image = bigImage(ROI.min(2):ROI.max(2), ROI.min(1):ROI.max(1), :);
[xnd,map] = rgb2ind(image,3,'nodither');
img3Colors = ind2rgb(xnd, map);
R = sum(sum(img3Colors(:,:,1)));
G = sum(sum(img3Colors(:,:,2)));
B = sum(sum(img3Colors(:,:,3)));
M = max([R,G,B]);
N = min([R,G,B]);
code = 1;
[corners, boardSize] = detectCheckerboardPoints(img3Colors, 'HighDistortion',true,'PartialDetections',true, 'MinCornerMetric',0.18);
if ((boardSize(1,1)) < 2 || (boardSize(1,2)<2))
    code = 0;
    squares=[];
    return;
end

x = reshape(corners(:,1),boardSize(1)-1,[]);
y = reshape(corners(:,2),boardSize(1)-1,[]);
if (max(x(1,:)) - min(x(1,:))) < (max(x(:,1)) - min(x(:,1)))
    x = x';
    y = y';
end
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
if (verbose)
    hold(ax,'on');
end
for row=1:NRows-1
    for col=1:NCols-1
        quad = [x(row,col)+1, y(row,col)+1; x(row+1, col)+1, y(row+1, col); x(row+1, col+1), y(row+1, col+1); x(row, col+1), y(row, col+1)+1];
        if (sum(any(isnan(quad))) == 0)
            [value, background] = GetSquareValues(image, imageBW, quad,prevailingColor, verbose);
            if (verbose)
                fprintf('col=%d row=%d color=%d value=%d, background=%d\n', col, row, prevailingColor, value, background);
            end
            quad(:,1) = quad(:,1) + ROI.min(1);
            quad(:,2) = quad(:,2) + ROI.min(2);
            if (ax ~= 0)
                scatter(ax, quad(:,1), quad(:,2),50,'oy', 'filled');
                drawnow;
            end
            squares(row,col).color =prevailingColor;
            squares(row,col).value =value;
            squares(row,col).quad = quad;
            squares(row,col).background = background;
        end
    end
end
blacks=[];
whites=[];
rotated = 0;
if (size(squares,1) > size(squares,2))
    % vertical array
    col = round(size(squares,2) / 2);
    for row=1:size(squares,1)
        if (squares(row,col).background == 1)
            whites = [whites; squares(row,col).value];
        end
    end
    rotated = (length(whites) > 1) && (sum(diff(whites)) == 0);
else
    % horizontal array
    row = round(size(squares,1) / 2);
    for col=1:size(squares,2)
        if (squares(row,col).background == 0)
            blacks = [blacks; squares(row,col).value];
        end
    end
    rotated = (length(blacks) > 1) && (sum(diff(blacks)) == 0);
end
%% Set Row value as per the whites
if (rotated)
    for col=1:NCols-1
        whites=[];
        for row=1:NRows-1
            if(squares(row, col).background == 1)
                whites = [whites; squares(row, col).value];
            end
        end
        if (size(whites,1)>0)
            rowNumber = mode (whites);
            for row=1:NRows-1
                squares(row, col).row = rowNumber;
                if (sum(whites) ~= length(whites)*rowNumber) || (rowNumber == 0)
                    squares(row, col).value = 0;
                end
    
            end
        end
    end
else
    for row=1:NRows-1
        whites=[];
        for col=1:NCols-1
            if(squares(row,col).background == 1)
                whites = [whites; squares(row,col).value];
            end
        end
        if (size(whites,1)>0)
            rowNumber = mode (whites);
            for col=1:NCols-1
                squares(row,col).row = rowNumber;
                if (sum(whites) ~= length(whites)*rowNumber) || (rowNumber == 0)
                    squares(row,col).value = 0;
                end
    
            end
        end
    end
end
%% Set Col value as per the blacks
if (rotated)
    for row=1:NRows-1
        blacks=[];
        for col=1:NCols-1
            if(squares(row,col).background == 0)
                blacks = [blacks; squares(row,col).value];
            end
        end
        if (size(blacks,1)>0)
            rowNumber = mode (blacks);
            for col=1:NCols-1
                squares(row,col).col = rowNumber;
                if (sum(blacks) ~= length(blacks)*rowNumber) || (rowNumber == 0)
                    squares(row,col).value = 0;
                end
            end
        end
    end
else
    for col=1:NCols-1
        blacks=[];
        for row=1:NRows-1
            if(squares(row,col).background == 0)
                blacks = [blacks; squares(row,col).value];
            end
        end
        if (size(blacks,1)>0)
            rowNumber = mode (blacks);
            for row=1:NRows-1
                squares(row,col).col = rowNumber;
                if (sum(blacks) ~= length(blacks)*rowNumber) || (rowNumber == 0)
                    squares(row,col).value = 0;
                end
            end
        end
    end
end
squares(:,end)=[];
squares(end,:)=[];
end


function [value, background] = GetSquareValues(image, imageBW, square, prevailingColor, verbose)
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
croppedBW = imageBW(y1:y2,x1:x2,:);
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
ci = ci(1:h,1:w,:);
cbw = imwarp(croppedBW,tf, OutputView = imref2d(size(imageBW))); 
cbw = cbw(1:h,1:w);
r=ci(:,:,1);
g=ci(:,:,2);
b=ci(:,:,3);
r(cbw)=0;
g(cbw)=0;
b(cbw)=0;
mo=cat(3,r,g,b);
background = sum(sum(sum(mo))) > 1500000;
% rad = roudnd((x2-x1)/16);
% se = strel('diamond',rad);
% sd = strel('diamond', rad-1);
% ei = imerode(cbw,se);
% di = imdilate(ei,sd);
dotAreas = regionprops('table', cbw,'Area');
da = (dotAreas.Area>7);
value = sum(da);
end

