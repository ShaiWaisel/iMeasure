close all; clear all;clc; imtool close all;
VERBOSE = 1;
platformID= "BaseP001";
configFileName = 'C:\Projects\iMeasure\FromErez\Database\!iMeasure.db';
dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SamCore1.db';
plyFileName = 'combined.ply';
SENSORS = 3;
FRAME = 1;

% Read sensors IDs and transformation matrices from configuration DB
connPlatform = sqlite(configFileName); 
platforms = sqlread(connPlatform,'Platforms');
sqlQuery = sprintf('SELECT * FROM Platforms WHERE PlatformID = "%s"', platformID);
platformData=table2array(fetch(connPlatform, sqlQuery));
close(connPlatform);

% Allocate memory
sensors = cell(1, SENSORS);
transMats = cell(1, SENSORS);
rawData = cell(1, SENSORS);
wall = cell(1, SENSORS);
transWall = cell(1, SENSORS);

% Get data from DB structure
for i=1:SENSORS
    sensors{i} = platformData(5+i*4);
    transMats{i} = sscanf(platformData(8+i*4), '%f;', [4, inf])';
end

% Open data DB and read raw data
conn = sqlite(dataFileName); 
for sensor=1:length(sensors)
    rawData{sensor} = sqlread(conn,sprintf('%s_Frames',sensors{sensor}));
    wall{sensor}= [];
    transWall{sensor}= [];
end
close(conn);

%create frame and points cloud from every frame
for sensor=1:SENSORS
    [wall{sensor}, u, v, img] = GrabFrameData(rawData{sensor}, FRAME, 0.0);
end

%transform raw points into platform coordinates
for sensor=1:SENSORS
    mat = transMats{sensor}; 
    tform = affinetform3d(mat');
    transWall{sensor} = pctransform(wall{sensor}, tform);
end

%show the result
player = pcplayer([-1.5 1.5], [-1.5 1.5], [0 1]);
walls = pccat([transWall{1} transWall{2} transWall{3}]);
view(player, walls);
pcwrite(walls,plyFileName);
return;


% ===================================================================================================

function [ptc, fu, fv, image] = GrabFrameData(data, frame, filterMargins)
    rgbW = 1280;
    rgbH = 720;
% build image raw data from DB 
    rgb = data(frame,19);
    rgb = cell2mat(rgb.FrameColor);
    rgb = reshape(rgb,3,[]);
    r = reshape(rgb(1,:), rgbW, rgbH)'; 
    g = reshape(rgb(2,:), rgbW, rgbH)'; 
    b = reshape(rgb(3,:), rgbW, rgbH)';
    image = cat(3, r, g, b);
    %figure, imshow(image);
%read XYZUV data from DB
    x = (cell2mat(data(frame, 20).FrameX));
    y = (cell2mat(data(frame, 21).FrameY));
    z = (cell2mat(data(frame, 22).FrameZ));
    u = (cell2mat(data(frame, 23).FrameU));
    v = (cell2mat(data(frame, 24).FrameV));
    fx=typecast(x,'single');
    fy=typecast(y,'single');
    fz=typecast(z,'single');
    fu=typecast(u,'single');
    fv=typecast(v,'single');
%filter out-of-image UV indices    
    idx = (fu<filterMargins) | (fu>(1.0-filterMargins)) | (fv<filterMargins) | (fv>(1.0-filterMargins));
    fx(idx)=[];
    fy(idx)=[];
    fz(idx) = [];
    fu(idx) = [];
    fv(idx) = [];
    iu = int32(floor(size(image,2)*fu)+1);
    iv = min(size(image,1),int32(floor(size(image,1)*fv)+1));
    linIdx = sub2ind(size(image), iv, iu);
%create points cloud in platform coordinates Z <= -Y, Y <= Z    
    ptc = pointCloud([fx, fz, -fy]);
    ptc.Color=[r(linIdx), g(linIdx), b(linIdx)];
end

