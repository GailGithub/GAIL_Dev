%This mfile may be used to download and install GAIL into 
%   the location you choose
%
%   Step 1.  Place this .m file where you want GAIL to go
%
%   Step 2.  Run this mfile

%% Download the package and change the directory
disp('The GAIL package is now downloading ...')
unzip('https://gail.googlecode.com/files/GAIL_1.0.0.zip') %download and unzip
cd('GAIL_1.0.0') %get to the right subdirectory
cd('GAIL_Matlab') %get to the right subdirectory

%% Install Gail
GAIL_Install %this installs GAIL
fprintf('\n\nNext we will run a quick test. Press return to continue...\n')
pause

%% Run a quick test
fprintf('\nmuhat=meanMC_g(@(n) rand(n,1))\n')
muhat=meanMC_g(@(n) rand(n,1)) %run meanMC_g once
fprintf('\nThis answer should be close to 0.5.\n\n')
%   Then you should be ready to use GAIL
disp('Next README.txt will be displayed.  Press return to continue...')
pause

%% Printing out README.txt
more on
type('README.txt')
more off