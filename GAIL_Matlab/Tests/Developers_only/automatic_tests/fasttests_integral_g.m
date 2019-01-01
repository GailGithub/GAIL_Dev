% fasttests_integral_g: fast tests for integral_g

%% CALL DOCTESTS
tic; doctest integral_g; time=toc
tic; doctest dt_integral_g ; time=toc

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
    warning('off','GAIL:integral_g:peaky')
    warning('off','GAIL:integral_g:exceedbudget')
    run_handle_ut('ut_integral_g')
    warning('on','GAIL:integral_g:peaky')
    warning('on','GAIL:integral_g:exceedbudget')
end
