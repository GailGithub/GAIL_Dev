% This m file is to download yuhanqeq package into the location you
% choose.
%
%   Step 1.  Place this file in the folder you want to use the package
%
%   Step 2.  Run this mfile
%% Download the package and change the directory
unzip('https://docs.google.com/uc?id=0B1VVw8pWbfwtREVkMWpfbXZvTmM&export=download') %download and unzip
cd('yuhanqeq') %get to the right subdirectory

qeqversion = 1;
fprintf('\nWelcome to yuhanqeq version %g.\n', qeqversion);
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
addpath(qeqp);% Add GAIL directories and subdirectories
savepath; % Save the changes
fprintf('\nyuhanqeq version %g has been installed successfully.\n\n', qeqversion);

reply = input('\nDo you want to run script to test yuhanqeq?\n(Your answer is not case sensitive.) Y/N [N]:','s');
if any(strcmpi(reply,{'yes','y'}));
workout_yuhanqeq
end

