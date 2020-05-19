function CopyZipGAILRepoAndGoogleDrive(newpath,GoogleDrivePath,zipfilename,excludefilelist)
% This script copies and combines files in both repository and GAIL's
% Google drive to make a release zip package.
%
% newpath --- path to where the zip is created, e.g., 'E:\GAIL_2_3\' for
% GAIL version 2.3
%
% GoogleDrivePath --- where your Google drive path, e.g.,
% 'E:\GoogleDrive\GAIL_Dev\'
%
% zipfilename --- the name of zip file, e.g., 'GAIL_2_3.zip' for GAIL 2.3
%
% excludefilelist --- name of an external file that contains the list of
% files we want to exclude in the zip, e.g., 'exclude2_3.lst' for GAIL 2.3
%
% CopyZipGAILRepoAndGoogleDrive('E:\GAIL_2_3\','E:\GoogleDrive\GAIL_Dev\','GAIL_2_3.zip','exclude2_3.lst')

GAILPATH = GAILstart(0); % Get the path of GAIL_Matlab
GoogleDrivePath1 = horzcat(GoogleDrivePath,'GAIL_Matlab',filesep);

%copy files in repository to new path
copyfile(GAILPATH,newpath,'f') % copy files in Google drive to new path
copyfile(GoogleDrivePath1,newpath,'f')

% zip files recursively excluding hidden file such as .ds_store and .icon
% in current folder
system(horzcat('zip -r ',zipfilename,' ',newpath,' -x@',excludefilelist));

movefile(zipfilename, GoogleDrivePath,'f') % move zip to Google drive path
display(horzcat(zipfilename,' is already put in ',GoogleDrivePath))
