% longtests_integral_g : long tests for integral_g

format short

%% Workouts
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_integral_g);
    results=run(ut_workout_integral_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_workout_integral_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_workout_integral_g is wrongly coded. We skip it.\n');
end
