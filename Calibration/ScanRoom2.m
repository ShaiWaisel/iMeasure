close all; clear all; imtool close all;
VERBOSE = 1;
platformID= "BaseP001";
%configFileName = 'C:\Projects\iMeasure\FromErez\Database\!iMeasure.db';
configFileName = 'C:\Projects\iMeasure\FromErez\Database\!iMeasure_LRB.db';
%dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SAMCore7114d.db';
%dataFileName = 'C:\Projects\iMeasure\FromErez\Database\SAMCore_frameframe_13_11_22_4.db';
dataFileName = 'C:\Projects\iMeasure\FromErez\Database\ScanLRB180_1.db';
plyFileName = 'combined.ply';

[sensors, transMats, rawData, platformData] = ReadFromSQL(configFileName, dataFileName, platformID);

figure; ax{1} = gca;
figure; ax{2} = gca;
figure; ax{3} = gca;
figure; ax{4} = gca;
scannedModel = StitchFrames(rawData, sensors, transMats,[1, size(rawData{1},1)], ax);
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

