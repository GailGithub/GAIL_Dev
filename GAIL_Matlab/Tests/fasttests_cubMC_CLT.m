% fasttests_cubMC_CLT: fast tests for cubMC_CLT

%% CALL DOCTESTS
tic; doctest cubMC_CLT; time = toc

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
    try
        Tests = matlab.unittest.TestSuite.fromClass(?e);
        results=run(ut_cubMC_CLT);
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        disp('Error: Test ut_cubMC_CLT is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_cubMC_CLT is wrongly coded. We skip it.\n');
    end    
    disp(results);
end