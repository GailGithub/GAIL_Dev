% fasttests_funmin_g_CSC
format
doctest funmin_g_CSC
doctest sinetest
doctest quadratic
doctest flat_bottom
doctest flatbottom2
doctest fejer_jackon_inequality

    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g_CSC);
        results=run(ut_funmin_g_CSC)
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        display('Error: Test ut_funmin_g is wrongly coded. We skip it.')
    end
    
    
 results=run(ut_workout_funmin_g_CSC)   