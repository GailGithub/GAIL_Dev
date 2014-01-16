%GAIL_Reinstall.m
% This script is to reinstall GAIL package. 
% To remove the old GAIL path and add the new one.
clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nYou are reinstalling GAIL version %g.\n\n',GAILVERSION);

if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
else
  gailp=genpath(GAILPATH); % Generate strings of paths to GAIL subdirectories 
end
rmpath(gailp);  % Remove GAIL paths from MATLAB search path.
addpath(gailp); % Add GAIL paths to MATLAB search path.
savepath;       % Save the changes.

fprintf('\nGAIL version %g has been reinstalled successfully.\n\n',GAILVERSION);