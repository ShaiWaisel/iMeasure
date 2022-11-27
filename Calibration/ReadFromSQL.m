
function [sensors, transMats, rawData, platformData] = ReadFromSQL(configFileName, dataFileName, platformID)
% Read sensors IDs and transformation matrices from configuration DB
connPlatform = sqlite(configFileName); 
platforms = sqlread(connPlatform,'Platforms');
sqlQuery = sprintf('SELECT * FROM Platforms WHERE PlatformID = "%s"', platformID);
platformData=table2array(fetch(connPlatform, sqlQuery));
close(connPlatform);
SENSORS = 0;
for i=1:4
    id = platformData(5+i*4);
    SENSORS = SENSORS + 1 - contains(id,"f000000");
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
for sensor=1:length(sensors)
    rawData{sensor} = sqlread(conn,sprintf('%s_Frames',sensors{sensor}));
    wall{sensor}= [];
    transWall{sensor}= [];
end
close(conn);
end

