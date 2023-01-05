function  SaveConfiguration(sensorsDBfile, platformID, transMats)
connPlatform = sqlite(sensorsDBfile); 
for camera=1:length(transMats)
    sqSensorName = sprintf('Sensor%d_tForm', camera);
    mat = transMats{camera}; 
    calibMats{camera} = mat;
    sqMat='';
    for i=1:4
        for j=1:4
            sqMat = strcat(sqMat, sprintf('%f; ',mat(i,j)));
        end
    end
    sqlOp =['UPDATE Platforms' , ...
        sprintf(' SET %s = "%s" ', sqSensorName, sqMat), ...
        sprintf(' WHERE PlatformID = "%s"', platformID)];
    exec(connPlatform, sqlOp);
end
close(connPlatform);
