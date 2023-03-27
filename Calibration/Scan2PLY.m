function Scan2PLY(varargin)
close all; 
clearvars -except varargin; 
clc;

if (length(varargin)< 4)
    disp(sprintf("Usage: Scan2PLY <scanFileName> <parametersFileName> <outputFileName> <verbose [0/1]\n>"));
    exit(-1);
end;
inputFileName = varargin(1);
disp(inputFileName)

if (~exist(inputFileName{1}, 'file'))
    disp(sprintf("Error: scan file [%s] not found!\n", inputFileName{1}));
    exit(-2);
end
parametersFileName = varargin(2);
disp(parametersFileName)
if (~exist(parametersFileName{1}, 'file'))
    disp(sprintf("Error: parameters file [%s] not found!\n", parametersFileName{1}));
    exit(-3);
end
outputFileName = varargin(3);
disp(outputFileName)
verbose = varargin(4);
dVerbose = str2double(verbose{1});
disp(sprintf("verbose: %d\n",dVerbose));
if (dVerbose>0)
    disp('starting app');
    system(sprintf("Scan2PLYApp %s %s %s%s",inputFileName{1},parametersFileName{1},outputFileName{1},"1"));
else
    disp('running silent')
    DoScan2PLY(inputFileName{1},parametersFileName{1},outputFileName{1}, 0);
end
exit(0);