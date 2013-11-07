% This mfile is to download and install SilvaQEQ package into the location 
% you choose.
%
%   Step 1.  Place this .m file where you want SilvaQEQ to go
%
%   Step 2.  Run this mfile
%% Download the package and change the directory
disp('The SilvaQEQ package is now downloading ...')
unzip('https://docs.google.com/uc?export=download&id=0ByjgFHKSc6HcUHF1UGUtTjZ1UTg') %download and unzip
cd('silvaQEQ') %get to the right subdirectory     

Silvaqeqversion = 1.0;
fprintf('\nWelcome to SilvaQeq version %g. \n', Silvaqeqversion);
Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    SilvaQEQPath=[pwd, PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    SilvaQEQPath = [pwd, PATHNAMESEPARATOR];
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
  qeqp=genpath(SilvaQEQPath);% adding all subdirectories 
end
addpath(qeqp);% Add GAIL directories and subdirectories
savepath; % Save the changes
fprintf('\nSilvaQeq version %g has been installed successfully.\n\n', Silvaqeqversion);

reply = input('\n Do you want to run script to test SilvaQeq?\n y/n [n]:','s');
if any(strcmpi(reply,{'yes','y'}));
workout_SilvaQEQ
end
