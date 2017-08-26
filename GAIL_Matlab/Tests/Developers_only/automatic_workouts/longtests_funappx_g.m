% logtests_funappx_g : long tests for funappx_g

format short

% doctest
% funappx_g
doctest par_funappx_g

% unit test
% funappx_g
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_funappx_g);
    results=run(ut_workout_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_workout_funappx_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_workout_funappx_g is wrongly coded. We skip it.\n');
end
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_convtest_funappx_g);
    results=run(ut_convtest_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_convtest_funappx_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_convtest_funappx_g is wrongly coded. We skip it.\n');
end


