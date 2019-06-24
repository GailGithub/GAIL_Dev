function fullfilename = save_mat(subdir, filename, varargin)
% SAVE_MAT: Save data to a MAT file in a subdirectory in 'OutputFiles'
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

if nargin <= 2
    error('Input values for isTimeStamped and variables to save.');
elseif nargin == 3 
    error('Input variables to save.');
else 
    isTimeStamped = varargin{1};
end

GAILPATH = GAILstart(0);
outputfolder =  [GAILPATH,'OutputFiles',filesep,subdir];
if exist(outputfolder,'dir') ~= 7
  mkdir(outputfolder);
end
if isTimeStamped == true
    fullfilename = strcat(outputfolder, filesep', filename,'-',...
        datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.mat');
else
    fullfilename = strcat(outputfolder, filesep', filename,'.mat');
end
varnames=cell(1,length(varargin)-1);
for k = 1:(nargin-3)
    varname = inputname(k+3);
    eval([varname, ' = varargin{k+1};']);
    varnames{k} = varname;
end
clear subdir filename GAILPATH varname outputfolder isTimeStamped k varargin;
if ~isempty(varnames)
   save(fullfilename, varnames{:});
else 
   warning('No variables saved.')
end
clear varnames;


