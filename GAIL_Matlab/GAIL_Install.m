% GAIL_INSTALL   Install GAIL. Add GAIL paths to MATLAB search path.
%clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nWelcome to GAIL version %g.\n\n', GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
end
gailp=genpath(GAILPATH);% Generate strings of paths to GAIL subdirectories
addpath(gailp);           % Add GAIL directories and subdirectories
savepath;                 % Save the changes
reply = input('\nDo you want to install html documentation files?\n (Your answer is not case sensitive.) Y/N [N]: ', 's');
if any(strcmpi(reply,{'yes','y'}));
    GAILpublish;        
    builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));
    fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
    addpath(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'helpsearch'));
else
    fprintf('\nYou can use test documentation to learn how to use GAIL.\n');
end
fprintf('\nGAIL version %g has been installed successfully.\n\n', GAILVERSION);

