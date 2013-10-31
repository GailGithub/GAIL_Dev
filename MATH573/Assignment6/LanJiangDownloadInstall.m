% This m file is to download LanJiang's qeq package into the location you
% choose.
%
%   Step 1.  Place this .m file in where you want LanJiangqeq to go
%
%   Step 2.  Run this mfile
%% Download the package and change the directory
unzip(['https://docs.google.com/uc?id=0B5as28avI00tdnE0UGFBMEx1NmM&export=download']) %download and unzip
cd('LanJiangqeqHW6') %get to the right subdirectory

LanJiangqeqversion = 1;
fprintf('\nWelcome to LanJiangqeq version %g.\n', LanJiangqeqversion);
Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    LANJIANGQEQPATH=[pwd, PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    LANJIANGQEQPATH = [pwd, PATHNAMESEPARATOR];
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
  qeqp=genpath(LANJIANGQEQPATH);% adding all subdirectories 
end
addpath(qeqp);% Add GAIL directories and subdirectories
savepath; % Save the changes
fprintf('\nLanJiangqeq version %g has been installed successfully.\n\n', LanJiangqeqversion);

reply = input('\nDo you want to run script to test LanJiangqeq?\n(Your answer is not case sensitive.) Y/N [N]:','s');
if any(strcmpi(reply,{'yes','y'}));
LanJiangTestqeq
end

