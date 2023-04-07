close all; clear all; imtool close all;
VERBOSE = 1;
platformID= "H001";
%configFileName = 'C:\Projects\iMeasure\FromErez\Database\!iMeasure.db';
configFileName = 'C:\Projects\iMeasure\FromErez\Database\#iMeasure.db';
%dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SAMCore7114d.db';
%dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SAMCore_frameframe_13_11_22_4.db';
dataFileName = 'C:\Projects\iMeasure\FromErez\Database\H001_Calib_22_11_22.db';
plyFileName = 'combined.ply';

[sensors, transMats, rawData, platformData] = ReadFromSQL(configFileName, dataFileName, platformID);

params.FrameGauge.Value = 0;
params.PlaneError = 0.02;
params.PlanePoints = 20000;
params.Margins = 0.1;
params.Axes{1} = gca;
params.Axes{2} = gca;
params.Axes{3} = gca;
params.Axes{4} = gca;
scannedModel = StitchFrames(rawData, sensors, transMats,[1, size(rawData{1},1)], params);
%show the result
player = pcplayer([-1.5 1.5], [-1.5 1.5], [-1 1]);
view(player.Axes, 10 , 80);

% Customize player axes labels
xlabel(player.Axes, 'X (m)')
ylabel(player.Axes, 'Y (m)')
zlabel(player.Axes, 'Z (m)')
hold(player.Axes, 'on');

view(player, scannedModel);
pcwrite(scannedModel,plyFileName);
return;

