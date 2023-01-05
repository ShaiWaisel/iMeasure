
function [code, sensors, transMats, rawData, platformData] = ReadFromSQL(configFileName, dataFileName, platformID)
% Read sensors IDs and transformation matrices from configuration DB
code = 0;
sensors = 0;
transMats = 0; 
rawData = 0;
platformData = 0;
if (~exist(configFileName, 'file'))
    code = -1;
    return;
end
if (~exist(dataFileName, 'file'))
    code = -2;
    return;
end
connPlatform = sqlite(configFileName); 
platforms = sqlread(connPlatform,'Platforms');
sqlQuery = sprintf('SELECT * FROM Platforms WHERE PlatformID = "%s"', platformID);
platformData=table2array(fetch(connPlatform, sqlQuery));
close(connPlatform);
SENSORS = 0;
for i=1:4
    id = platformData(5+i*4);
    SENSORS = SENSORS + 1 - contains(id,"f000000") - (id == "<empty>");
end

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
try
    for sensor=1:length(sensors)
        rawData{sensor} = sqlread(conn,sprintf('%s_Frames',sensors{sensor}));
        wall{sensor}= [];
        transWall{sensor}= [];
    end
catch
    code = -4;
end
close(conn);
end

