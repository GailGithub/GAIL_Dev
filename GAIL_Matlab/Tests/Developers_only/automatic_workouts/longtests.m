% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,MATLABVERSION]  = GAILstart(false);
ver_str = sprintf('MatlabVer%2.2f', MATLABVERSION);
ver_str = strrep(ver_str,'.','-')
filename = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_workouts-', datestr(now,'yyyy-mm-dd-HH-MM-SS'), '_', ver_str, '.txt');
diary(filename)
tic;

format short

% print matlab version
ver

% print chebfun version
ver_chebfun = ver('Chebfun');
if ~isempty(ver_chebfun)
  ver_chebfun
else
  ver_matlab = ver;
  index = find(strcmp({A.Name}, 'Chebfun')==1);
  if index > 0 
      fprintf('\nChebfun version is %s.\n\n', ver_matlab(index).Version);
  else
      fprintf('\nChebfun version not found.\n');
  end
end

%% utilties
run_handle_ut('ut_save_mat')
if usejava('jvm') || MATLABVERSION <= 7.12
  run_handle_ut('ut_save_eps')
end

%% Workouts
fprintf('=============================Starting longtests at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_cubLattice_g
fprintf('=============================Finished longtests_cubLattice_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_cubMC_g
fprintf('=============================Finished longtests_cubMC_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_cubSobol_g
fprintf('=============================Finished longtests_cubSobol_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_funappx_g
fprintf('=============================Finished longtests_funappx_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_funmin_g
fprintf('=============================Finished longtests_funmin_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_integral_g
fprintf('=============================Finished longtests_integral_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
longtests_meanMC_g
fprintf('=============================Finished longtests_meanMC_g at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))
%% Papers
%longtests_conePaper

%% doctests and unit tests for deprecated algos
longtests_deprecated_algos

fprintf('=============================Finished all the longtests at %s \n', datestr(now,'yyyy-mm-dd-HH-MM-SS'))

time=toc;
disp(time)

diary off
