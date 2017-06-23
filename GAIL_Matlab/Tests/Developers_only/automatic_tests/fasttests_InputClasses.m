% fasttests_InputClasses: fast tests for input classes

%% CALL DOCTESTS
tic; doctest gail.gail_in_param; time=toc
tic; doctest gail.gail1D_in_param; time=toc
tic; doctest gail.gailMD_in_param; time=toc

%% CALL UNIT TESTS