% fasttests: Drives all fast doctests and unit tests
format short

completereport = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_tests-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

[GAILPATH,GAILVERSION,MATLABVERSION] = GAILstart(false)
%shortutestreport = strcat(GAILPATH,'OutputFiles',filesep,...
%    'gail_unittests','.txt');
%fid = fopen(shortutestreport,'wt');

tstart=tic;

%% CALL fasttests for individual algorithms
fasttests_integral_g
fasttests_funappx_g
fasttests_funmin_g
fasttests_meanMC_g
fasttests_cubMC_g
fasttests_cubLattice_g
fasttests_cubSobol_g
fasttests_meanMC_CLT

%% CALL fasttests for other key GAIL components
format short

fasttests_OptionPricing
fasttests_InputClasses

total_time=toc(tstart)
 
diary off
%fclose(fid);
format
