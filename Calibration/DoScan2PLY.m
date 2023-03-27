function DoScan2PLY(scanFile, paramsFile, outputFile, verbose)
[code, ~, sensors, ~,  transMats, rawData] = ReadFromSQLscan(scanFile);
if (code ==0)
    fid=fopen(paramsFile);
    C = textscan(fid,'%f');
    fclose(fid);
    C=C{1};
    params.fromFrame = C(1);
    params.toFrame = C(2);
    params.PlaneError = C(3);
    params.ICPgrid=C(4);
    params.PointsTolerance = C(5);
    params.Margins = C(6);
    params.PlanePoints = C(7);
    params.AlignModel = C(8);
    frameSpan = [C(1), C(2)];
    scannedModel = StitchFrames(rawData, sensors, transMats, frameSpan, params, verbose);
    pcwrite(scannedModel,outputFile);
end
