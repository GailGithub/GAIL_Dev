function save_eps(subdir, filename)
% SAVE_MAT: Save figure to an eps file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   subdir          string, name of subdirectory 
%   filename        string, filename of eps file
%   
% Example:
%   save_eps('ConesPaperOutput','ConesPaperFunAppxTest', ...
%   'tauvec','pini','pfin','succnowarn', 'succwarn','failnowarn','failwarn');

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
path = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    subdir, PATHNAMESEPARATOR, filename, '-',...
    datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.eps');
print('-depsc', path)