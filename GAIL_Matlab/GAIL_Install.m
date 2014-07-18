% GAIL_INSTALL   Install GAIL. Add GAIL paths to MATLAB search path.
clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nWelcome to GAIL version %g.\n', GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
end
gailp=genpath(GAILPATH);% Generate strings of paths to GAIL subdirectories
cd(strcat(GAILPATH,'Utilities'));
GAILpublish; % Publish searchable html documentation
docPath = strcat(GAILPATH,'Documentation');
cd(docPath);
helpPath = strcat(docPath,PATHNAMESEPARATOR,'html');
builddocsearchdb(helpPath);
addpath(gailp);           % Add GAIL directories and subdirectories
savepath;                 % Save the changes
fprintf('\nGAIL version %g has been installed successfully.\n\n', GAILVERSION);