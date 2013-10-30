% UT_JAGADEESQEQ unit test for jagadeeQeq
classdef ut_jagadeesQeq < matlab.unittest.TestCase
    
    methods(Test)
        % ut_jagadeesQeq tests solutions to the quadratic equation
        % a*xˆ2 + b*x + c = 0

        function QeqOrdered(testCase)
            tol = 1e-4;
            solNow = jagadeesQeq(1,2,3);
            exactSol = [(-1.0000 -1.4142i),  (-1.0000 + 1.4142i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function QeqVector(testCase)
            tol = 1e-4;
            solNow = jagadeesQeq([9,5,3]);
            exactSol = [(-0.2778 - 0.5061i),  (-0.2778 + 0.5061i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end

        function QeqLabel(testCase)
            tol = 1e-4;
            solNow = jagadeesQeq('a', 8, 'c', 45, 'b', 2);
            exactSol = [(-0.1250 - 2.3684i), ( -0.1250 + 2.3684i)];
            error = norm(exactSol-solNow);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
end
% 
% run(ut_jagadeesQeq)
% Running ut_jagadeesQeq
% ...
% Done ut_jagadeesQeq
% __________
% 
% 
% ans = 
% 
%   1x3 TestResult array with properties:
% 
%     Name
%     Passed
%     Failed
%     Incomplete
%     Duration
% 
% Totals:
%    3 Passed, 0 Failed, 0 Incomplete.
%    0.011467 seconds testing time.

