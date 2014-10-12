% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
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
warning('off','MATLAB:integral_g:exceedbudget')
warning('off','MATLAB:integral_g:peaky')
run_handle('workout_integral_g')
warning('on','MATLAB:integral_g:peaky')
warning('on','MATLAB:integral_g:exceedbudget')


%% Papers
% Cone paper
run_handle('ConesPaperFoolFunctions')

% MCQMC paper
run_handle('RunTestcubMConGeoAsianCall')
run_handle('RunTestcubMConGaussian')
run_handle('RunTestcubMConGaussiand1')

% Function minimization
run_handle('workout_testerrortolerance')
run_handle('workout_testXtolerance')
run_handle('workout_testboth')
run_handle('workout_twoextreme')

run_handle('UniFunMin_Plot_Bump')
run_handle('UniFunMin_Plot_Flat')
run_handle('UniFunMin_Plot_TwoExtreme')
run_handle('UniFunMin_test_ErrorTolerance')
run_handle('UniFunMin_test_ErrorXTolerance')
run_handle('UniFunMin_test_TwoExtreme')
run_handle('UniFunMin_test_XTolerance')
 

%% Unit tests
MATLABVERSION = gail.matlab_version
if MATLABVERSION >= 8  
    Tests = matlab.unittest.TestSuite.fromClass(?ut_ConesPaper);
    run(ut_ConesPaper)
end

time=toc;
disp(time)

diary off
