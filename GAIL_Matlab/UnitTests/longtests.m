% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_workouts-', datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.txt');
diary(filename)
tic; 

%% Workouts 
% doctests
format short
doctest dt_meanMC_g_TrafficModel
Test_MeanMC_g
doctest dt_cubMC_g
test_cubMC_g

% other workouts
warning('off','MATLAB:integraltau_g:peaky')
run_handle('tryout_integral_g_FJH')
warning('on','MATLAB:integraltau_g:peaky')
warning('off','MATLAB:integral_g:exceedbudget')
warning('off','MATLAB:integral_g:peaky')
run_handle('workout_integral_g')
warning('on','MATLAB:integral_g:peaky')
warning('on','MATLAB:integral_g:exceedbudget')


%% Papers
% Cone paper
run_handle('ConesPaperFoolFunctions')
run_handle('conepaper_test_integral_g')
run_handle('conepaper_test_funappx_g')

% MCQMC paper
run_handle('RunTestcubMConGeoAsianCall')
run_handle('RunTestcubMConGaussian')
run_handle('RunTestcubMConGaussiand1')

% Function minimization
run_handle('workout_testerrortolerance')
run_handle('workout_testXtolerance')
run_handle('workout_testboth')
run_handle('workout_twoextreme')

%% Unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8  
    %run(ut_funappx_g)   
end

time=toc;
disp(time)

diary off
