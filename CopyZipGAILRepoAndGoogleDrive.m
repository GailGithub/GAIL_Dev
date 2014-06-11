% This script is to copy and combine files in both repository and Google drive
% to make a release zip package.
clear all; close all; clc;
newfolder = 'GAIL_1_3';
newpath = horzcat('E:\',newfolder);
% set up a new path where you want to copy things
GAILPATH = GAILstart;%Get the path of repository from GAILstart
GoogleDrivePath1 = 'E:\GoogleDrive\GAIL_Dev\GAIL_Matlab';
GoogleDrivePath = 'E:\GoogleDrive\GAIL_Dev';
% get the google drive path.
copyfile(GAILPATH,newpath)%copy files in repository to new path
copyfile(GoogleDrivePath1,newpath)% copy files in google drive to new path
% zip files recursively excluding hidden file such as .ds_store and .icon
% in current folder
system(horzcat('zip -r GAIL_1_3.zip ', newpath, ' -x@exclude1_3.lst'));  
% move zip to google drive path
movefile('GAIL_1_3.zip', GoogleDrivePath,'f')