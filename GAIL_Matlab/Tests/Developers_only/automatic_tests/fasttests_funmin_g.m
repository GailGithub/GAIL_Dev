% fasttests_funmin_g: fast tests for funmin_g

%% CALL DOCTESTS
tic; doctest funmin_g; time=toc


%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
  run_handle_ut('ut_funmin_g')
  run_handle_ut('ut_funmin_g_end')
end