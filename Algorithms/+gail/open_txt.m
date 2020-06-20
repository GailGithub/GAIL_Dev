function [fileID, fullfilename] = open_txt(subdir, filename, isTimeStamped)
% OPEN_TXT: Opens a text file in a subdirectory in 'OutputFiles'
% 
% Inputs:
%   subdir          name of subdirectory 
%   filename        filename of mat file
%   isTimeStamped   boolean, set to true to include a timestamp in filename
%   variable_names  variables in workspace of calling function to persist
%   
% Example:
%   save_mat('ConesPaperOutput','ConesPaperFunAppxTest', ...
%      tauvec,pini,pfin,succnowarn);

if nargin < 3
    isTimeStamped = true;
elseif nargin < 2 
    error('Not enough inputs.');
end

GAILPATH = GAILstart(0);
outputfolder =  [GAILPATH,'OutputFiles',filesep,subdir];
if exist(outputfolder,'dir') ~= 7,
  mkdir(outputfolder);
end
if isTimeStamped == true,
    fullfilename = strcat(outputfolder, filesep', filename,'-',...
        datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
else
    fullfilename = strcat(outputfolder, filesep', filename,'.txt');
end
fileID = fopen(fullfilename,'w');



