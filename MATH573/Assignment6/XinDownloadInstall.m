% This m file will download and install Xinqeq package into
% the location you choose.
%
% Step 1. Place this folder where you want
%
% Step 2. Run this m.file

%% Download the package and change the directory
disp('The Xinqeq package is now downloading ...')
unzip('https://docs.google.com/uc?id=0B7DrUtLuOcgkS245TVQ3X0RIVk0&export=download') %download and unzip
cd('Xinqeq') %get to the right subdirectory

Xinqeqversion=1;
fprintf('\nWelcome to Xinqeq version %g.\n', Xinqeqversion);
Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    XINQEQPATH = [pwd, PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    XINQEQPATH = [pwd, PATHNAMESEPARATOR];
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
  Xinqeqp=genpath(XINQEQPATH);% adding all subdirectories 
end
addpath(Xinqeqp);
savepath; % Save the changes
fprintf('\nXinqeq version %g has been installed successfully.\n\n', Xinqeqversion);



%% Run a quick test
fprintf('The solution of quadratic equation  x^2 + 4x + 3 = 0')
x=Xinqeq(1,4,3) % run Xinqeq once
fprintf('\nThis answer should be x =  -3  -1 \n\n') % Then you should be ready to use GAIL


