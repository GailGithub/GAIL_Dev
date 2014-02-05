% Drives all lengthy doctests and unit tests

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_workouts-', datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.txt');
diary(filename)

tic; 
% Call doctest 
format short
doctest dt_meanMC_g_TrafficModel
doctest dt_cubMC_g

Test_MeanMC_g
test_cubMC_g
warning('off','MATLAB:integraltau_g:peaky')
tryout_integral_g_FJH
warning('on','MATLAB:integraltau_g:peaky')
warning('off','MATLAB:integral_g:exceedbudget')
warning('off','MATLAB:integral_g:peaky')
workout_integral_g 
warning('on','MATLAB:integral_g:peaky')
warning('on','MATLAB:integral_g:exceedbudget')

% Call unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8  
    %run(ut_funappx_g)   
end

time=toc;

disp(time)

diary off
