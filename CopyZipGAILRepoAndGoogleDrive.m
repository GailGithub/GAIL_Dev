% This script is to copy and combine files in both repository and Google drive
% to make a release zip package.
clear all; close all; clc;
newpath = '/Users/BabyAudrey/Documents/Study/Research/GAIL_1_3_0';
% set up a new path where you want to copy things
GAILPATH = GAILstart;%Get the path of repository from GAILstart
GoogleDrivePath = '/Users/BabyAudrey/Google Drive/GAIL_Dev/GAIL_Matlab';
% get the google drive path.
copyfile(GAILPATH,newpath)%copy files in repository to new path
copyfile(GoogleDrivePath,newpath)% copy files in google drive to new path
% zip files recursively excluding hidden file such as .ds_store and .icon 
cd('/Users/BabyAudrey/Documents/Study/Research/')
system('zip -r GAIL130.zip GAIL_1_3_0 -x@exclude1_3.lst');  %-x@exclude1_3.lst
%zip -r GAIL130.zip GAIL_1_3_0 -x@exclude1_3.lst
%system(horzcat('zip -r GAIL_1_3_0 ', newpath, ' -x@exclude1_3.lst'));  
