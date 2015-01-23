function fullfilename = save_mat(subdir, filename,varargin)
% SAVE_MAT: Save data to a MAT file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   subdir          name of subdirectory 
%   filename        filename of mat file
%   variable_names  names of variables in work space to save
%   
% Example:
%   save_mat('ConesPaperOutput','ConesPaperFunAppxTest', ...
%   'tauvec','pini','pfin','succnowarn', 'succwarn','failnowarn','failwarn');

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
path = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    subdir, PATHNAMESEPARATOR');
fn = strcat(filename,'-',datestr(now,'yyyy-mmm-dd-HH-MM-SS'),'.mat');
% fullfilename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
%     subdir, PATHNAMESEPARATOR', filename,'-',...
%     datestr(now,'yyyy-mmm-dd-HH-MM-SS'),'.mat');
fullfilename = [path,fn];
save([path,fn],varargin)
%, varargin)