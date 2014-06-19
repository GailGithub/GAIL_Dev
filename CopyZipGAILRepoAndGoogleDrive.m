% This script is to copy and combine files in both repository and Google drive
% to make a release zip package.
function CopyZipGAILRepoAndGoogleDrive(newpath,GoogleDrivePath,zipfilename,excludefilelist)
% newpath is where you want to copy things including folder name
% eg. 'E:\GAIL_Version#'
% GoogleDrivePath is where your google drive path
% eg. 'E:\GoogleDrive\GAIL_Dev'
% Get the path of repository from GAILstart
% zipfilename is the name of zip file
% eg. 'GAIL_Version#'
% excludefilelist is the list of file we want to exclude
% eg. 'excludeVersion#.lst'
[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart;
% Get the path where GAIL_Matlab locate
GoogleDrivePath1 = horzcat(GoogleDrivePath,PATHNAMESEPARATOR,'GAIL_Matlab');
%copy files in repository to new path
copyfile(GAILPATH,newpath)
% copy files in google drive to new path
copyfile(GoogleDrivePath1,newpath)
% zip files recursively excluding hidden file such as .ds_store and .icon
% in current folder
system(horzcat('zip -r ',zipfilename,' ',newpath, ' -x@',excludefilelist));  
% move zip to google drive path
movefile(zipfilename, GoogleDrivePath,'f')
display(horzcat(zipfilename,' is already put in ', GoogleDrivePath))