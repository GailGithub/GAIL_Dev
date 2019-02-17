% fasttests_funappx_g: fast tests for funappx_g

%% CALL DOCTESTS 
tic; doctest funappx_g; time=toc
tic; doctest dt_funappx_g; time=toc

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
  run_handle_ut('ut_funappx_g')
end