% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION]  = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_workouts-', datestr(now,'yyyymmdd'),'.txt');
diary(filename)
tic; 

format short

%% Workouts 
% meanMC_g
doctest dt_meanMC_g_TrafficModel
run_handle('Test_meanMC_g')

% cubMC_g
doctest dt_cubMC_g
run_handle('Test_cubMC_g')

% integral_g
run_handle('workout_integral_g')


%% Papers
% Cone paper
if usejava('jvm')
  run_handle('ConesPaperFoolFunctions')
end
if MATLABVERSION >= 8  
    run(ut_ConesPaper)
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
if MATLABVERSION >= 8  
    run(ut_UniFunMin)
end

time=toc;
disp(time)

diary off
