% UT_YUHANQEQ unit test for yuhanqeq
classdef ut_yuhanqeq < matlab.unittest.TestCase
    
    methods(Test)
        
        function yuhanqeqof231(testCase)
            tol = 1e-6;
            yuhansol = yuhanqeq([2,3,1]);
            exactsol= [-1 -0.5] ;
            actualerr = norm(exactsol-yuhansol);
            testCase.verifyLessThanOrEqual(actualerr,tol);
        end
        
        function yuhanqeqof121(testCase)
            yuhansol = yuhanqeq(1,-2,1);
            exactsol= [1 1] ;
            testCase.verifyEqual(yuhansol,exactsol);
        end
        
        function yuhanqeqof012(testCase)
            yuhansol = testCase.verifyWarning(@()yuhanqeq(0,1,2),'MATLAB:yuhanqeq:zeroquadterm');
            exactsol= -2;
            testCase.verifyEqual(yuhansol,exactsol);
        end
        
    end
end

% Running ut_yuhanqeq
% ...
% Done ut_yuhanqeq
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
%    2.1596 seconds testing time.
