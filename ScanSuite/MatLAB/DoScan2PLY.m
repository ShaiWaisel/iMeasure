function [scannedModel] = DoScan2PLY(scanFile, paramsFile, outputFile, Axes, FrameGauge, DisplayGauge, verbose)
[code, ~, sensors, ~,  transMats, rawData] = ReadFromSQLscan(scanFile);
scannedModel=[];
if (code ==0)
    params = ReadScan2PLYparams(paramsFile);
    params.Axes = Axes;
    params.FrameGauge = FrameGauge;
    params.DisplayGauge = DisplayGauge;
    frameSpan = [params.fromFrame, params.toFrame];
    scannedModel = StitchFrames(rawData, sensors, transMats, frameSpan, params, verbose);
    pcwrite(scannedModel,outputFile);
end
end