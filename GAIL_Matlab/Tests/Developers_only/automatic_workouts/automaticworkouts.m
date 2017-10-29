%% Set the directory for running our matlab test
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/
% Install the latest
GAIL_Install

%% In Karlin server, chebfun is installed at this location
% cd /home/gail/GAIL_tests
% git clone https://github.com/chebfun/chebfun.git
chebfunroot = '/home/gail/GAIL_tests/chebfun';
addpath(chebfunroot), savepath

% Go to the tests direcoty (% is comment in matlab environment)
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_workouts/

fprintf('Launching longtests at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))

longtests

fprintf('Exiting longtests at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))

exit
