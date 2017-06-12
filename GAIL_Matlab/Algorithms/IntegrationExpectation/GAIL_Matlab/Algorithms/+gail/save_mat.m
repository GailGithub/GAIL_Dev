function fullfilename = save_mat(subdir, filename, varargin)
% SAVE_MAT: Save data to a MAT file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   subdir          name of subdirectory 
%   filename        filename of mat file
%   variable_names  variables in workspace to save
%   
% Example:
%   save_mat('ConesPaperOutput','ConesPaperFunAppxTest', ...
%   tauvec,pini,pfin,succnowarn');

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
outputfolder =  [GAILPATH,'OutputFiles',PATHNAMESEPARATOR,subdir];
if exist(outputfolder) ~= 7,
  mkdir(outputfolder);
end
fullfilename = strcat(outputfolder, PATHNAMESEPARATOR', filename,'-',...
    datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.mat');
varnames={};
for k = 1:length(varargin)
    varname = inputname(k+2);
    eval([varname, ' = varargin{k};']);
    varnames{k} = varname;
end
clear subdir filename GAILPATH PATHNAMESEPARATOR varname outputfolder
save(fullfilename, varnames{:});