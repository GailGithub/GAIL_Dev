% fasttests_funmin_g: fast tests for funmin_g

%% CALL DOCTESTS
tic; doctest funmin_g; time=toc


%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g);
        results=run(ut_funmin_g)
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        display('Error: Test ut_funmin_g is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_funmin_g is wrongly coded. We skip it.\n');
    end
    
    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g_end);
        results=run(ut_funmin_g_end)
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        display('Error: Test ut_funmin_g_end is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_funmin_g_end is wrongly coded. We skip it.\n');
    end
end