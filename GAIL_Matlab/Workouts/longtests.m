% Drives all lengthy doctests and unit tests

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_workouts-', datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.txt');
diary(filename)
tic; 
%% Workouts 

warning('off','MATLAB:integraltau_g:peaky')
tryout_integral_g_FJH
warning('on','MATLAB:integraltau_g:peaky')
warning('off','MATLAB:integral_g:exceedbudget')
warning('off','MATLAB:integral_g:peaky')
workout_integral_g 
warning('on','MATLAB:integral_g:peaky')
warning('on','MATLAB:integral_g:exceedbudget')

format short
doctest dt_meanMC_g_TrafficModel
Test_MeanMC_g
doctest dt_cubMC_g
test_cubMC_g


%% Papers
ConesPaperFoolFunctions
conepaper_test_integral_g					
conepaper_test_funappx_g



%% Unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8  
    %run(ut_funappx_g)   
end

time=toc;

disp(time)

diary off
