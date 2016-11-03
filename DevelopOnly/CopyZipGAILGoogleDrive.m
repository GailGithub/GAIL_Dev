function CopyZipGAILGoogleDrive(GoogleDrivePath,zipFile)
% This script copies files in GAIL's Google drive to make a release zip package.
%
% GoogleDrivePath --- where your Google drive path, e.g.,
% 'E:\GoogleDrive\GAIL_Dev\'
%
% zipFile --- the name of zip file, e.g., 'GAIL_2_0.zip' for GAIL 2.0
%
% CopyZipGAILGoogleDrive('E:\GoogleDrive\GAIL_Dev','GAIL_Build.zip')

currentPath = pwd;
cd(GoogleDrivePath);
system(horzcat('zip -r ',zipFile,' GAIL_Matlab -x@excludeGoogle.m'));
cd(currentPath);
display(horzcat(zipFile,' is already put in ',GoogleDrivePath));
end