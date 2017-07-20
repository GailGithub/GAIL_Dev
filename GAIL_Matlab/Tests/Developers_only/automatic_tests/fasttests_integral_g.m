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
    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_integral_g);
        results=run(ut_integral_g)
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        display('Test ut_integral_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_integral_g is wrongly coded. We skip it.\n');
    end
    warning('on','GAIL:integral_g:peaky')
    warning('on','GAIL:integral_g:exceedbudget')
end
