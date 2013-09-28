% GAIL_Install.m
% This script is to install GAIL, to add GAIL path to MATLAB search path.
clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nWelcome to GAIL version %g.\n', GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
else
  gailp=genpath(GAILPATH);% adding all subdirectories 
end
addpath(gailp);% Add GAIL directories and subdirectories
savepath; % Save the changes
fprintf('\nGAIL version %g has been installed successfully.\n\n', GAILVERSION);