% This script is to copy and combine files in both repository and Google drive
% to make a release zip package.
clear all; close all; clc;
newpath = '/Users/BabyAudrey/Documents/Study/Research/GAIL_1_3_0';
% set up a new path where you want to copy things
GAILPATH = GAILstart;
%Get the path of repository from GAILstart
GoogleDrivePath = '/Users/BabyAudrey/Google Drive/GAIL_Dev/GAIL_Matlab';
% get the google drive path.
copyfile(GAILPATH,newpath)%copy files in repository to new path
copyfile(GoogleDrivePath,newpath)% copy files in google drive to new path
delete('/Users/BabyAudrey/Documents/Study/Research/GAIL_1_3_0/.*')
zip('GAIL_1_3_0',newpath)% zip it.
%zip -r GAIL_1_3_0.zip bitvolution -x *.DS_Store*
% Lan is still working on how to excluded hidden file such as .dsstore and
% .icon in every single folder. please give her some time.
