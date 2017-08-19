% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,MATLABVERSION]  = GAILstart(false);
filename = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_workouts-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(filename)
tic;

format short

%% utilties
run_handle_ut('ut_save_mat')
if usejava('jvm') || MATLABVERSION <= 7.12
  run_handle_ut('ut_save_eps')
end

%% Workouts
longtests_cubLattice_g
longtests_cubMC_g
longtests_cubSobol_g
longtests_funappx_g
longtests_funmin_g
longtests_integral_g
longtests_meanMC_g

%% Papers
longtests_conePaper

%% doctests and unit tests for deprecated algos
longtests_deprecated_algos



time=toc;
disp(time)

diary off
