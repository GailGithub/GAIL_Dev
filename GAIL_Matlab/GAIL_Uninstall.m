% GAIL_Uninstall.m
% This script is to uninstall GAIL. 
% Remove path from Matlab search path and delete all the files and folders.
clear all;close all;clc;
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart(0);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
else
  gailp=genpath(GAILPATH);
end
fprintf('\nYou are going to remove GAIL v %g path.\n', GAILVERSION);
reply = input('\nDo you want to remove GAIL path from MATLAB search path?\n(Your answer is not case sensitive.) Y/N [N]:','s');
if any(strcmpi(reply,{'yes','y'}));
    rmpath(gailp);% Remove path from MATLAB search path
    savepath; % Save the changes
    fprintf('\nGAIL path has been removed from MATLAB search path successfully.\n')
    fprintf('\nYou are going to delete GAIL v %g files.\n ', GAILVERSION);
    reply = input('\nDo you want to delete ALL GAIL files?\n (Your answer is not case sensitive.) Y/N [N]: ', 's');
    if any(strcmpi(reply,{'yes','y'}));
        delete('GAILPATH/*.*') % delete all GAIL files
        rmdir(GAILPATH,'s')
        fprintf('\nAll GAIL v %g files have been deleted successfully.\n', GAILVERSION);
    else
        fprintf('\nGAIL files have not been deleted.\n');
    end
else
    fprintf('\nGAIL path has not been removed.\n');
end