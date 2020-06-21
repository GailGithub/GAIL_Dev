function DownloadInstallGail_2_3_1()
% This function downloads and installs GAIL version 2.3.1 into the
% location you choose
%
%   Step 1.  Place this M-file where you want GAIL to go
%
%   Step 2.  Run this M-file in MATLAB
%
%This file installs GAIL version 2.3.1 in the subdirectory
%"GAIL_2.3.1/GAIL_Matlab".
%% Download the package and change the directory
disp('The GAIL package is now being downloaded...')
unzip('http://www.math.iit.edu/~Meshfree-methods-seminar/GAIL/GAIL_2_3_1.zip', 'GAIL_2_3_1') %download and unzip
cd('GAIL_2_3_1')  %get to the right subdirectory
addpath(genpath(pwd));
%% Install GAIL
GAIL_Install      %this installs GAIL
fprintf('\n\nNext we will run a quick test. Press RETURN to continue...\n')
pause

%% Run a quick test
fprintf('\nmuhat=meanMC_g(@(n) rand(n,1))\n');
muhat=meanMC_g(@(n) rand(n,1)); %run meanMC_g once
fprintf('\nmuhat = %.4f\n',muhat);
fprintf('\nThis answer should be close to 0.5.\n\n');
%Then you should be ready to use GAIL
disp('Next README.txt will be displayed. Press RETURN to continue to ');
disp('the next line, or press the spacebar to advance to the next page...');
pause

%% Print out README.md
moreStatus = get(0,'More');
more on
type('README.md')
more(moreStatus);
end
