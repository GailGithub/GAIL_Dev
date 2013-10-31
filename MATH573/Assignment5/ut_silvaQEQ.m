% ut_silvaQEQ unit test for silvaQEQ
% Name: Tiago Silva
% Email: tsilva@hawk.iit.edu

classdef ut_silvaQEQ < matlab.unittest.TestCase
    
    methods(Test)
        
        function SilvaOrderedCoeffCase(testCase)
            tol = 1e-5;
            testSol = silvaQEQ(1,2,1);
            exactSol = [-1,  -1];
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
        function SilvaVectorCase(testCase)
            tol = 1e-5;
            testSol = silvaQEQ([1,-2,1]);
            exactSol = [1,  1];
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
        function SilvaVararginCase(testCase)
            tol = 1e-5;
            coeff.a = 1; coeff.b = 3; coeff.c = 2;
            testSol = silvaQEQ(coeff);
            exactSol = [-2, -1];
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
    
end


% Running ut_silvaQEQ
% ...
% Done ut_silvaQEQ
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
%    0.013428 seconds testing time.