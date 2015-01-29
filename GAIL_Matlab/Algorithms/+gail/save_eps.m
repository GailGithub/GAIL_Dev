function path = save_eps(subdir, filename)
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
if exist(subdir) ~= 7,
  mkdir(strcat([GAILPATH,'OutputFiles',PATHNAMESEPARATOR], subdir));
end
path = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    subdir, PATHNAMESEPARATOR, filename, '-',...
    datestr(now,'yyyy-mmm-dd-HH-MM-SS'),'.eps');
print('-depsc', path)