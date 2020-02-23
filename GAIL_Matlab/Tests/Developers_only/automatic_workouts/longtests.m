% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,MATLABVERSION]  = GAILstart(false);
ver_str = sprintf('MatlabVer%2.2f', MATLABVERSION);
ver_str = strrep(ver_str,'.','-')
filename = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_workouts-', datestr(now,'yyyy-mm-dd-HH-MM-SS'), '_', ver_str, '.txt');
diary(filename)
tstart=tic;

format short
format compact

% print gail version
disp(['GAILVERSION = ', num2str(GAILVERSION)]);

% print matlab version
ver

% print chebfun version
chebfun_ver = gail.ver('Chebfun')

%% utilties
run_handle_ut('ut_save_mat')
if usejava('jvm') || MATLABVERSION <= 7.12
    run_handle_ut('ut_save_eps')
end


%% Workouts
if license('test', 'Signal_Toolbox') && MATLABVERSION >= 9.7
    longtests_cubBayesNet_g
    fprintf('=============================Finished longtests_cubBayesNet_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
end
longtests_cubBayesLattice_g
fprintf('=============================Finished longtests_cubBayesLattice_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))

fprintf('=============================Starting longtests at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_cubLattice_g
fprintf('=============================Finished longtests_cubLattice_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_cubSobol_g
fprintf('=============================Finished longtests_cubSobol_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_meanMC_g
fprintf('=============================Finished longtests_meanMC_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_cubMC_g
fprintf('=============================Finished longtests_cubMC_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_funappx_g
fprintf('=============================Finished longtests_funappx_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_funmin_g
fprintf('=============================Finished longtests_funmin_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_integral_g
fprintf('=============================Finished longtests_integral_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_meanMC_CLT
fprintf('=============================Finished longtests_meanMC_CLT at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))

%% doctests and unit tests for deprecated algos
%longtests_deprecated_algos

fprintf('=============================Finished all the longtests at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))

total_time = toc(tstart)

diary off
