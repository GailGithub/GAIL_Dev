% fasttests_meanMC_CLT: fast tests for meanMC_CLT

%% CALL DOCTESTS
tic; doctest meanMC_CLT; time = toc

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMC_CLT);
        results=run(ut_meanMC_CLT);
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        disp('Error: Test ut_meanMC_CLT is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_meanMC_CLT is wrongly coded. We skip it.\n');
    end    
    disp(results);
end