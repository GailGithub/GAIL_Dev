function path = save_image(figh, subdir, filename, varargin)
% SAVE_IMAGE: Save figure as an image file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   figh            handle, handle to the figure
%   subdir          string, name of subdirectory 
%   filename        string, filename of eps file
%   isTimeStamped   boolean, default or set to false to include 
%                   a timestamp in the filename
% 
% Example:
%   save_image(figh,'ConesPaperOutput','ConesPaperFunAppxTest',true);

if nargin <= 3
    isTimeStamped = false;
else 
    isTimeStamped = varargin{1};
end
filetype = '.png';

GAILPATH = GAILstart(0);
dirpath=strcat([GAILPATH,'OutputFiles',filesep], subdir);
if exist(dirpath) ~= 7,
  mkdir(dirpath);
end
if isTimeStamped == true,
    path = strcat(dirpath, filesep, filename, '-',...
        datestr(now,'yyyy-mm-dd-HH-MM-SS'),filetype);
else
    path = strcat(dirpath, filesep, filename,filetype);
end

saveas(figh, path)
