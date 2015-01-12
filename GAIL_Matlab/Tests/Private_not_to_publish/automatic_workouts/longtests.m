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

% meanMCBer_g
run_handle('Test_MeanMCBer_g')

% cubMC_g
doctest dt_cubMC_g
run_handle('Test_cubMC_g')

% integral_g
run_handle('workout_integral_g')

%Workout of function approximation
if MATLABVERSION >= 8  
    run(ut_workoutfunappx_g)
end

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

% meanMCBer_g paper
if MATLABVERSION >= 8  
    run(ut_meanMCBer_g)
end

% cubSobol_g paper
run_handle('RunTestCubatureonGeoAsianCall')
run_handle('RunTestCubatureonKeister')

% cubLattice_g paper
run_handle('TestFourierTransform');
run_handle('test_cubLattice');
run_handle('RunTestCubatureonGeoAsianCall');
run_handle('RunTestCubatureonKeister')

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
