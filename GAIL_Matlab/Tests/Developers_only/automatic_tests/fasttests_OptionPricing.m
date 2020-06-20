% fasttests_OptionPricing: fast tests for OptionPricing objects

%% CALL DOCTESTS
tic; doctest assetPath; time=toc
tic; doctest optPayoff; time=toc
tic; doctest optPrice; time=toc
tic; doctest brownianMotion; time=toc
tic; doctest stochProcess; time=toc
tic; doctest whiteNoise; time=toc

%% CALL UNIT TESTS
run(ut_brownianMotion)
run(ut_stochProcess)
run(ut_whiteNoise)