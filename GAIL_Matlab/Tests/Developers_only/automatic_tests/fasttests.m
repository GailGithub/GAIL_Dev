% fasttests: Drives all fast doctests and unit tests
format short 
format compact

[GAILPATH,GAILVERSION,MATLABVERSION]  = GAILstart(false);
ver_str = sprintf('MatlabVer%2.2f', MATLABVERSION);
ver_str = strrep(ver_str,'.','-');
completereport = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_tests-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),  '_', ver_str,'.txt');
diary(completereport)
disp(['GAILVERSION = ', num2str(GAILVERSION),'.  MATLABVERSION = ', num2str(MATLABVERSION)]);


tstart=tic;

%% CALL fasttests for individual algorithms
if license('test', 'Signal_Toolbox') && MATLABVERSION >= 9.7
    fasttests_cubBayesNet_g
end
fasttests_cubBayesLattice_g
fasttests_integral_g
fasttests_funappx_g
fasttests_funmin_g
fasttests_meanMC_g
fasttests_cubMC_g
fasttests_cubLattice_g
fasttests_cubSobol_g


%% CALL fasttests for other key GAIL components
format short

fasttests_OptionPricing
fasttests_InputClasses

total_time = toc(tstart)
 
diary off
format
