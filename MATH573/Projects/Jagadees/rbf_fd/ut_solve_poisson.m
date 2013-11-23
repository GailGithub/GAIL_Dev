% UT_solve_poisson unit test for solve_poisson
classdef ut_solve_poisson < matlab.unittest.TestCase
    
    methods(Test)
        % ut_solve_poisson tests solutions to the RBF-FD Approximation

        function test1(testCase)
            tol = 1e-4;
            solNow = solve_poisson(1,2,3);
            exactSol = [(-1.0000 -1.4142i),  (-1.0000 + 1.4142i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function test2(testCase)
            tol = 1e-4;
            solNow = solve_poisson([9,5,3]);
            exactSol = [(-0.2778 - 0.5061i),  (-0.2778 + 0.5061i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function test3(testCase)
            tol = 1e-4;
            solNow = solve_poisson('a', 8, 'c', 45, 'b', 2);
            exactSol = [(-0.1250 - 2.3684i), ( -0.1250 + 2.3684i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
end
