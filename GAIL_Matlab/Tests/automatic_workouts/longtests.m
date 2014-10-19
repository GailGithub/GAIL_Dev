% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION]  = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_workouts-', datestr(now,'yyyymmdd'),'.txt');
diary(filename)
tic; 

%% Workouts 
% doctests
format short
doctest dt_meanMC_g_TrafficModel
Test_meanMC_g
doctest dt_cubMC_g
Test_cubMC_g

% other workouts
run_handle('workout_integral_g')


%% Papers
% Cone paper
if usejava('jvm')
  run_handle('ConesPaperFoolFunctions')
end

% MCQMC paper
run_handle('RunTestcubMConGeoAsianCall')
run_handle('RunTestcubMConGaussian')
run_handle('RunTestcubMConGaussiand1')

% Function minimization
if usejava('jvm')
  run_handle('UniFunMin_Plot_Bump')
  run_handle('UniFunMin_Plot_Flat')
  run_handle('UniFunMin_Plot_TwoExtreme')
end
 

%% Unit tests
if MATLABVERSION >= 8  
    run(ut_ConesPaper)
    run(ut_UniFunMin)
end

time=toc;
disp(time)

diary off
