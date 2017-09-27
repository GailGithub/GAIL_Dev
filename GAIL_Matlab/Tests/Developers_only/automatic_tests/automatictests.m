%% Set the directory for running our matlab test
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/
% Install the latest
GAIL_Install
% Go to the tests direcoty (% is comment in matlab environment)
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/
%runtests %> /home/gail/GAIL_tests/test_results.txt former tests
fasttests
exit
