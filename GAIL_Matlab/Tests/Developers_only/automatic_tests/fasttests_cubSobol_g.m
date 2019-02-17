% fasttests_cubcubSobol_g: fast tests for cubSobol_g

%% CALL DOCTESTS
tic; doctest cubSobol_g; time=toc


%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
  run_handle_ut('ut_cubSobol_g')
end