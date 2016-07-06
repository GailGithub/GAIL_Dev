function GAIL_Install(isBuild,isPublish)
%GAIL_INSTALL Install or reinstall GAIL. Add GAIL paths to MATLAB search path.
%
% isBuild: true to download and explode a zip file. Default to false.
% isPublish: true to publish HTML format of GAIL documentation. Default to false.
%
if nargin < 1
    isBuild =  false;
end
if nargin < 2
    isPublish = false;
end
[GAILPATH,GAILVERSION,MATLABVERSION] = GAILstart(false);
fprintf('\nWelcome to Guaranteed Automatic Integration Library (GAIL).\nYou are installing GAIL v%s.\n\n', GAILVERSION);
if MATLABVERSION < 7,
    error('This version is only supported on MATLAB 7 and above.');
end
gailp=genpath(GAILPATH);% Generate strings of paths to GAIL subdirectories
warninfo = warning('query','MATLAB:rmpath:DirNotFound');
warning('off',warninfo.identifier);
rmpath(gailp);% Remove path from MATLAB search path
warning(warninfo.state, warninfo.identifier);
addpath(gailp);           % Add GAIL directories and subdirectories
savepath;                 % Save the changes
if isBuild
   unzip('http://www.math.iit.edu/~Meshfree-methods-seminar/GAIL/GAIL_Build.zip',[GAILPATH,'..',filesep]);
end
if isPublish
   gail.GAIL_Publish(isPublish,~isPublish,isPublish);
end
fprintf('\nGAIL version %s has been installed successfully.\n\n', GAILVERSION);
end