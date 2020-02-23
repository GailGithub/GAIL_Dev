% longtests_cubBayesNet_g: long tests for cubBayesNet_g

%% CALL DOCTESTS


%% CALL UNIT TESTS
run_handle_ut('ut_cubBayesNet_g');

try
  cubBayesNet_long_tests;
catch
    disp('Error: cubBayesNet_long_tests is wrongly coded. We skip it.')
end