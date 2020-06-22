function path = save_eps(subdir, filename, varargin)
% SAVE_EPS: Save figure to an eps file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   subdir          string, name of subdirectory 
%   filename        string, filename of eps file
%   isTimeStamped   boolean, default or set to true to include 
%                   a timestamp in the filename
% 
% Example:
%   save_eps('ConesPaperOutput','ConesPaperFunAppxTest', ...
%   'tauvec','pini','pfin','succnowarn', 'succwarn','failnowarn','failwarn');

if nargin <= 2 
    isTimeStamped = true;
else 
    isTimeStamped = varargin{1};
end

GAILPATH = GAILstart(0);
dirpath=strcat([GAILPATH,'OutputFiles',filesep], subdir);
if exist(dirpath) ~= 7,
  mkdir(dirpath);
end
if isTimeStamped == true,
    path = strcat(dirpath, filesep, filename, '-',...
        datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.eps');
else
    path = strcat(dirpath, filesep, filename,'.eps');
end

print('-depsc', path)