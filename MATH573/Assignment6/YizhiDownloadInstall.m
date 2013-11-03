%This mfile may be used to download and install YizhiQEQ into 
%   the location you choose
%
%   Step 1.  Place this .m file in where you want YizhiQEQ to go
%
%   Step 2.  Run this mfile

%% Download the package and change the directory
unzip('https://sites.google.com/site/yizhiqeq/home/downloads/YizhiQEQ.zip')
disp('Download is in process...')
cd('YizhiQEQ') %get to the right subdirectory


%% Install YizhiQEQ
YIZHIQEQVERSION = 1.0;

Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    YIZHIQEQPATH=[pwd, PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    YIZHIQEQPATH = [pwd, PATHNAMESEPARATOR];
else
    error('I don''t recognize this computer.')
end

V = version;
wh = strfind(V,'.');
wh = wh(2)-1;
MATLABVERSION = str2double(V(1:wh));

fprintf('\nWelcome to YizhiQEQ version %g.\n', YIZHIQEQVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
else
  yizhiqeqp=genpath(YIZHIQEQPATH);% adding all subdirectories 
end
addpath(yizhiqeqp);% Add YizhiQEQ directories and subdirectories
savepath; % Save the changes
fprintf('\nYIZHIQEQ version %g has been installed successfully.\n\n', YIZHIQEQVERSION);

%% Run a quick test
disp('Running a test...')
result = YizhiQEQ([1,-4,3])
disp('The results above should be to 1 and 3.')
%   Then you should be ready to use YizhiQEQ