
function [code, sensors, devices, locations, transMats, rawData] = ReadFromSQLscan(dataFileName)
% Read sensors IDs and transformation matrices from configuration DB
code = 0;
devices = 0;
sensors = 0;
transMats = 0; 
locations = 0;
rawData = 0;
if (~exist(dataFileName, 'file'))
    code = -1;
    return;
end
conn = sqlite(dataFileName); 
platforms = sqlread(conn, 'DEVICEMAPPING');

SENSORS = size(platforms,1);
% Allocate memory
sensors = cell(1, SENSORS);
devices = cell(1, SENSORS);
locations = cell(1, SENSORS);
transMats = cell(1, SENSORS);
rawData = cell(1, SENSORS);
wall = cell(1, SENSORS);
transWall = cell(1, SENSORS);

% Get data from DB structure
for i=1:SENSORS
    sensors{i} = platforms{i,1};
    devices{i} = platforms{i,2};
    locations{i} = platforms{i,4};
    transMats{i} = eye(4);
end

% Open data DB and read raw data
try
    for device=1:length(devices)
        rawData{device} = sqlread(conn,sprintf('%s_Frames',devices{device}));
        wall{device}= [];
        transWall{device}= [];
    end
catch
    code = -2;
end
close(conn);
end

