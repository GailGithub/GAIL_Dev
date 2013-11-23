% UT_stencil_support_selection unit test for stencil_support_selection
classdef ut_stencil_support_selection < matlab.unittest.TestCase
    
    methods(Test)
        % ut_stencil_support_selection tests solutions to the stencil_support_selection

        function test1(testCase)
            tol = 1e-4;
            solNow = stencil_support_selection(1,2,3);
            exactSol = [(-1.0000 -1.4142i),  (-1.0000 + 1.4142i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function test2(testCase)
            tol = 1e-4;
            solNow = stencil_support_selection([9,5,3]);
            exactSol = [(-0.2778 - 0.5061i),  (-0.2778 + 0.5061i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function test2(testCase)
            tol = 1e-4;
            solNow = stencil_support_selection('a', 8, 'c', 45, 'b', 2);
            exactSol = [(-0.1250 - 2.3684i), ( -0.1250 + 2.3684i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
end

