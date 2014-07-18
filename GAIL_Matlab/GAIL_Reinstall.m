%GAIL_REINSTALL  Reinstall GAIL. Remove existing GAIL paths and add new ones.
clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nYou are reinstalling GAIL version %g.\n\n',GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
end
gailp=genpath(GAILPATH); % Generate strings of paths to GAIL subdirectories 
rmpath(gailp);  % Remove GAIL paths from MATLAB search path.
addpath(gailp); % Add GAIL paths to MATLAB search path.
savepath;       % Save the changes.
cd(strcat(GAILPATH,'Utilities'));
GAILpublish; % Publish searchable html documentation
docPath = strcat(GAILPATH,'Documentation');
cd(docPath);
helpPath = strcat(docPath,PATHNAMESEPARATOR,'html');
builddocsearchdb(helpPath);
cd('..');
fprintf('\nGAIL version %g has been reinstalled successfully.\n\n',GAILVERSION);