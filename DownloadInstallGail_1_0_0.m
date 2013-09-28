%This mfile may be used to download and install GAIL into 
%   the location you choose
%
%   Step 1.  Place this .m file in where you want GAIL to go
%
%   Step 2.  Run this mfile

%% Download the package and change the directory
unzip(['http://math.iit.edu/~openscholar/sites/default/files/meshfree/'...
    'files/gail_1.0.0.zip']) %download and unzip
cd('GAIL_1.0.0') %get to the right subdirectory
cd('GAIL_Matlab') %get to the right subdirectory

%% Install Gail
GAIL_Install %this installs GAIL

%% Run a quick test
muhat=meanMC_g(@(n) rand(n,1)); %run meanMC_g once
disp(['The number, ' num2str(muhat,2) ', should be close to 0.5'])
%   Then you should be ready to use GAIL