% UT_find_stencil_weights unit test for find_stencil_weights
classdef ut_find_stencil_weights < matlab.unittest.TestCase
    
    methods(Test)
        % ut_find_stencil_weights tests find_stencil_weights

        function test1(testCase)
            tol = 1e-4;
            solNow = find_stencil_weights(1,2,3);
            exactSol = [(-1.0000 -1.4142i),  (-1.0000 + 1.4142i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function test2(testCase)
            tol = 1e-4;
            solNow = find_stencil_weights([9,5,3]);
            exactSol = [(-0.2778 - 0.5061i),  (-0.2778 + 0.5061i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function test2(testCase)
            tol = 1e-4;
            solNow = find_stencil_weights('a', 8, 'c', 45, 'b', 2);
            exactSol = [(-0.1250 - 2.3684i), ( -0.1250 + 2.3684i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
end
