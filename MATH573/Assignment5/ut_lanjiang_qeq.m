classdef ut_lanjiang_qeq < matlab.unittest.TestCase
    
    methods(Test)
        
        function LanJiangqeqof23(testCase)
            tol = 1e-2;
            Lansol = LanJiangqeq([1,-5,6]);
            exactsol= [2 3] ;
            actualerr = norm(exactsol-Lansol);
            testCase.verifyLessThanOrEqual(actualerr,tol);
        end
        
        function LanJiangqeqof11(testCase)
            Lansol = LanJiangqeq(1,2,1);
            exactsol= [-1 -1] ;        
            testCase.verifyEqual(Lansol,exactsol);
        end
    end
end

% Running ut_lanjiang_qeq
% ..
% Done ut_lanjiang_qeq
% __________
% 
% 
% ans = 
% 
%   1x2 TestResult array with properties:
% 
%     Name
%     Passed
%     Failed
%     Incomplete
%     Duration
% 
% Totals:
%    2 Passed, 0 Failed, 0 Incomplete.
%    0.0044535 seconds testing time.