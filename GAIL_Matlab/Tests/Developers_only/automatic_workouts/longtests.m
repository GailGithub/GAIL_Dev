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
fprintf('Starting longtests at %s \n', datetime)
longtests_cubLattice_g
fprintf('Finished longtests_cubLattice_g at %s \n', datetime)
longtests_cubMC_g
fprintf('Finished longtests_cubMC_g at %s \n', datetime)
longtests_cubSobol_g
fprintf('Finished longtests_cubSobol_g at %s \n', datetime)
longtests_funappx_g
fprintf('Finished longtests_funappx_g at %s \n', datetime)
longtests_funmin_g
fprintf('Finished longtests_funmin_g at %s \n', datetime)
longtests_integral_g
fprintf('Finished longtests_integral_g at %s \n', datetime)
longtests_meanMC_g
fprintf('Finished longtests_meanMC_g at %s \n', datetime)
%% Papers
%longtests_conePaper

%% doctests and unit tests for deprecated algos
longtests_deprecated_algos

fprintf('Finished all the longtests at %s \n', datetime)

time=toc;
disp(time)

diary off
