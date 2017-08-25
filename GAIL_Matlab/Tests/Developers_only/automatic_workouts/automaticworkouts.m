%% Set the directory for running our matlab test
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/
% Install the latest
GAIL_Install

%% In Karlin server, chebfun is installed at this location
% cd /home/gail/GAIL_tests
% git clone https://github.com/chebfun/chebfun.git
chebfunroot = '/home/gail/GAIL_tests/chebfun';
addpath(chebfunroot), savepath
# print the version
help chebfun

% Go to the tests direcoty (% is comment in matlab environment)
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_workouts/

longtests
exit
