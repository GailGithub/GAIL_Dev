%GAIL_REINSTALL  Reinstall GAIL. Remove existing GAIL paths and add new ones.
%clear all; close all; clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nYou are reinstalling GAIL version %s.\n\n',GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
end
gailp=genpath(GAILPATH); % Generate strings of paths to GAIL subdirectories 
warninfo = warning('query','MATLAB:rmpath:DirNotFound');
warning('off',warninfo.identifier);
rmpath(gailp);% Remove path from MATLAB search path
warning(warninfo.state, warninfo.identifier);
addpath(gailp); % Add GAIL paths to MATLAB search path.
savepath;       % Save the changes.
reply = input('\nDo you want to install html documentation files?\n (Your answer is not case sensitive.) Y/N [N]: ', 's');
    if any(strcmpi(reply,{'yes','y'}));
        GAILpublish;        
        builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));
        fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
    else
        fprintf('\nYou can use test documentation to learn how to use GAIL.\n');
    end
fprintf('\nGAIL version %s has been reinstalled successfully.\n\n', GAILVERSION);
