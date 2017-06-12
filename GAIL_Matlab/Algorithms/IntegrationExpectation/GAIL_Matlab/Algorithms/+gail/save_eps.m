function path = save_eps(subdir, filename)
% SAVE_EPS: Save figure to an eps file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   subdir          string, name of subdirectory 
%   filename        string, filename of eps file
%   
% Example:
%   save_eps('ConesPaperOutput','ConesPaperFunAppxTest', ...
%   'tauvec','pini','pfin','succnowarn', 'succwarn','failnowarn','failwarn');

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
dirpath=strcat([GAILPATH,'OutputFiles',PATHNAMESEPARATOR], subdir);
if exist(dirpath) ~= 7,
  mkdir(dirpath);
end
path = strcat(dirpath, PATHNAMESEPARATOR, filename, '-',...
    datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.eps');
print('-depsc', path)