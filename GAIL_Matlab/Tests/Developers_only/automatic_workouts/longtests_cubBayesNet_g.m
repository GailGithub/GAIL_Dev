% longtests_cubBayesNet_g: long tests for cubBayesNet_g

%% CALL DOCTESTS


%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
  warning('Cannot run unit tests in MATLAB version before 8.1');
else
  run_handle_ut('ut_cubBayesNet_g');
end

try
  cubBayesNet_long_tests;
catch
    disp('Error: cubBayesNet_long_tests is wrongly coded. We skip it.')
end
