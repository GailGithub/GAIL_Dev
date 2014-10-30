%GAIL_REINSTALL  Reinstall GAIL. Remove existing GAIL paths and add new ones.
%clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nYou are reinstalling GAIL version %g.\n\n',GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
end
gailp=genpath(GAILPATH); % Generate strings of paths to GAIL subdirectories 
warninfo = warning('query','MATLAB:rmpath:DirNotFound');
warning('off',warinfo.identifier);
rmpath(gailp);% Remove path from MATLAB search path
warning(warninfo.state, warninfo.identifier);
addpath(gailp); % Add GAIL paths to MATLAB search path.
savepath;       % Save the changes.
warninfo = warning('query','MATLAB:doc:DocNotInstalled');
warning('off', warninfo.identifier);
builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));
warning(warninfo.state, warninfo.identifier);
fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
fprintf('\nGAIL version %g has been reinstalled successfully.\n\n', GAILVERSION);
