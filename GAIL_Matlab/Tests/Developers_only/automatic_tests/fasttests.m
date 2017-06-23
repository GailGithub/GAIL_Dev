% Drives all doctests and unit tests
format short

GAILPATH = GAILstart(0);
completereport = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_tests-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

shortutestreport = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_unittests','.txt');
fid = fopen(shortutestreport,'wt');

tic

%% CALL fasttests for individual algorithms
fasttests_integral_g
fasttests_funappx_g
fasttests_funmin_g
fasttests_meanMC_g
fasttests_cubMC_g
fasttests_cubLattice_g
fasttests_cubSobol_g

fasttests_OptionPricing


%tic; doctest gail.gail_in_param; time=toc
%tic; doctest gail.gail1D_in_param; time=toc
%tic; doctest gail.gailMD_in_param; time=toc


time=toc
 
diary off
fclose(fid);
format
