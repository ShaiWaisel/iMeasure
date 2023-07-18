function WriteScan2PLYparams(params, paramsFileName)
    fid=fopen(paramsFileName,'w');
    fprintf(fid,'%d\n',round(params.fromFrame));
    fprintf(fid,'%d\n',round(params.toFrame));
    fprintf(fid,'%f\n',(params.PlaneError));
    fprintf(fid,'%f\n',(params.ICPgrid));
    fprintf(fid,'%f\n',(params.PointsTolerance));
    fprintf(fid,'%f\n',params.Margins);
    fprintf(fid,'%d\n',round(params.PlanePoints));
    fprintf(fid,'%d\n',round(params.AlignModel));
    fclose(fid);
end
