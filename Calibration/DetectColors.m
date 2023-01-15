function [code, ROI, mask] = DetectColors(image, colorCode, minPixSport, verbose)
    R=image(:,:,1);
    G=image(:,:,2);
    B=image(:,:,3);
    g=rgb2gray(image);
    r1=imsubtract(R,g);
    g1=imsubtract(G,g);
    b1=imsubtract(B,g);
    code = 0;
    y1 = (r1>5)  & (g1 > 5);
    rm = max(max(r1));
    imsubtract(r1, uint8(y1).*rm);
    rg = max(max(g1));
    imsubtract(g1, uint8(y1).*rg);
    r1 = (r1>30);
    g1 = (g1>20);
    b1 = (b1>20);
    r1 = bwareaopen(r1, minPixSport);
    g1 = bwareaopen(g1, minPixSport);
    b1 = bwareaopen(b1, minPixSport);
    croppedImage = image;
    mask = im2bw(image);
    ROI.min = [0, 0];
    ROI.max = [0, 0];
    if (colorCode == 1) && (sum(sum(r1)) > 10000)
        code = 1;
        [ROI, mask] = Crop(image, r1);
    end

    if (colorCode == 2) && (sum(sum(g1)) > 10000)
        code = 1
        [ROI, mask] = Crop(image, g1);
    end

    if (colorCode == 3) && (sum(sum(b1)) > 10000)
        code = 1;
        [ROI, mask] = Crop(image, b1);
    end

    if (colorCode == 4) && (sum(sum(y1)) > 10000)
        code = 1;
        [ROI, mask] = Crop(image, y1);
    end

end

function [ROI, croppedMask] = Crop(image, mask)
p  = regionprops (mask, 'BoundingBox');
a = struct2array(p);
a = reshape(a,4,[])';
MIN=min(a(:,1:2));
r=[a(:,1)+a(:,3), a(:,2)+a(:,4)];
MAX=max(r);
ROI.min = [ceil(MIN(1)), ceil(MIN(2))];
ROI.max = [floor(MAX(1)), floor(MAX(2))];
croppedMask = mask(ROI.min(2):ROI.max(2),ROI.min(1):ROI.max(1), :);
%figure; imshow(r1); hold on;rectangle('Position',[MIN,MAX-MIN],'EdgeColor','y');

end