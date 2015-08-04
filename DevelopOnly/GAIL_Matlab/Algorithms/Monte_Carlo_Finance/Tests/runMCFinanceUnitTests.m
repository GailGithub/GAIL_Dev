function runMCFinanceUnitTests()

%% runMCFinanceUnitTests 
% runs unit tests of the programs in the Monte Carlo in Finance directory
% to make sure that all is working well.

runtests(fileparts(which('ut_stochProcess')))

