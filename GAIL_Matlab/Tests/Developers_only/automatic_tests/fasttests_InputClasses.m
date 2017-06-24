% fasttests_InputClasses: fast tests for input classes

%% CALL DOCTESTS
format short
tic; doctest gail.gail_in_param; time=toc
tic; doctest gail.gail1D_in_param; time=toc
tic; doctest gail.gailMD_in_param; time=toc
format long
tic; doctest gail.errorParam; time = toc
tic; doctest gail.fParam; time = toc
tic; doctest gail.cubParam; time = toc


%% CALL UNIT TESTS

