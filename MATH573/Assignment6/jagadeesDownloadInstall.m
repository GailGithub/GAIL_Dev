%This mfile may be used to download and install jagadeesQeq into
% the location you choose
%
% Step 1. Place this .m file where you want jagadeesQeq to go
%
% Step 2. Run this mfile
%% Download the package and change the directory
disp('The jagadeesQeq package is now downloading ...')
unzip('https://github.com/jagadeesr/M573/archive/jagadeesQeq_V1_0.zip') %download ...and unzip
cd('M573-jagadeesQeq_V1_0') %get to the right subdirectory

qeqversion = 1.0;
fprintf('\nWelcome to jagadeesQeq version %2.1f \n', qeqversion);
Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    QEQPATH=[pwd, PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    QEQPATH = [pwd, PATHNAMESEPARATOR];
else
    error('I don''t recognize this computer.')
end


V = version;
wh = strfind(V,'.');
wh = wh(2)-1;
MATLABVERSION = str2double(V(1:wh));
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
else
  qeqp=genpath(QEQPATH);% adding all subdirectories 
end
addpath(qeqp);% Add directories and subdirectories
savepath; % Save the changes
fprintf('\n jagadeesQeq version %g has been installed successfully.\n\n', qeqversion);


reply = input('\n Do you want to run script to test jagadeesQeq?\n y/n [n]:','s');
if any(strcmpi(reply,{'yes','y'}));
workout_jagadeesQeq
end

