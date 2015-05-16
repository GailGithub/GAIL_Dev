function GAIL_Install(isBuild)
%GAIL_INSTALL Install or reinstall GAIL. Add GAIL paths to MATLAB search path.
if nargin < 1
    isBuild =  false;
end
[GAILPATH,GAILVERSION,MATLABVERSION] = GAILstart;
fprintf('\nWelcome to Guaranteed Automatic Integration Library (GAIL).\nYou are installing GAIL v%s.\n\n', GAILVERSION);
if MATLABVERSION < 7,
    error('This version is only supported on Matlab 7.x and above.');
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
ifGenerateHtml=true; ifGenerateLateX=true; ifBuildSearchIndex=true;
gail.GAIL_Publish(ifGenerateHtml, ifGenerateLateX, ifBuildSearchIndex);
fprintf('\nGAIL version %s has been installed successfully.\n\n', GAILVERSION);
end