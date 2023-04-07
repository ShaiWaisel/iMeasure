function WriteScan2PLYparams(params, paramsFileName)
    fid=fopen(paramsFileName,'w');
    fprintf('%d\n',round(params.fromFrame));
    fprintf('%d\n',round(params.toFrame));
    fprintf('%f\n',(params.PlaneError));
    fprintf('%f\n',(params.ICPgrid));
    fprintf('%f\n',(params.PointsTolerance));
    fprintf('%d\n',round(params.Margins*100));
    fprintf('%d\n',round(params.PlanePoints));
    ffprintf('%d\n',round(params.AlignModel));
    fclose(fid);
end
