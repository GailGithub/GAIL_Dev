% longtests_cubBayesLattice_g: long tests for cubBayesLattice_g

%% CALL DOCTESTS
try
  doctest dt_cubBayesLattice_g
catch ME
    disp('Exception: longtests_cubBayesLattice_g : in doctest dt_cubBayesLattice_g. not ok.')
    msgText = getReport(ME); display(msgText)
end

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
  warning('Cannot run unit tests in MATLAB version before 8.1');
else
  run_handle_ut('ut_cubBayesLattice_g');
end

try
  cubBayesLattice_long_tests
catch err
  disp(['File: ', err.stack(1).name, '; Line number: ', int2str(err.stack(1).line),  '; Error: ',  err.message ])
  disp('Error: cubBayesLattice_long_tests is wrongly coded. not ok. We skip it.')
end
