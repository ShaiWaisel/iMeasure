function [scannedModel] = DoScan2PLY(scanFile, paramsFile, outputFile, Axes, FrameGauge, DisplayGauge, verbose)
[code, ~, sensors, ~,  transMats, rawData] = ReadFromSQLscan(scanFile);
scannedModel=[];
if (code ==0)
    if (size(rawData{1},1) > 0)
        params = ReadScan2PLYparams(paramsFile);
        if (params.fromFrame == 0)
            params.fromFrame = 1;
        end
        if (params.toFrame == 0)
            params.toFrame = size(rawData{1},1);
        end
        params.Axes = Axes;
        params.FrameGauge = FrameGauge;
        params.DisplayGauge = DisplayGauge;
        frameSpan = [params.fromFrame, params.toFrame];
        scannedModel = StitchFrames(rawData, sensors, transMats, frameSpan, params, verbose);
        pcwrite(scannedModel,outputFile);
    end
end
    end