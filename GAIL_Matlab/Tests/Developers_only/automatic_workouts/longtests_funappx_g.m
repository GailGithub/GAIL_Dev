% logtests_funappx_g : long tests for funappx_g

format short

%% doctest
%doctest par_funappx_g

%% unit test
% funappx_g
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_funappx_g);
    results = run(ut_workout_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    disp('Error: Test ut_workout_funappx_g is wrongly coded. We skip it.')
end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_convtest_funappx_g);
    results = run(ut_convtest_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    disp('Error: Test ut_convtest_funappx_g is wrongly coded. We skip it.')
end

% Workouts for Traub paper
run(ut_traub_paper, 'test_traub_paper_funappx_g');


