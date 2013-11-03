%UT_INTEGRAL_G unit test for integral_g
classdef ut_yizhiqeq < matlab.unittest.TestCase    
    methods (Test)
        
        function tauchange(testCase)
            coeff=[1,4,3];
            estroots=yizhiqeq(coeff);
            actroots=[-3,-1];
            testCase.verifyEqual(estroots,actroots);
        end
              
    end
end

%%Results
% >> run(ut_yizhiqeq)
% Running ut_yizhiqeq
% .
% Done ut_yizhiqeq
% __________
% 
% ans = 
%   TestResult with properties:
% 
%           Name: 'ut_yizhiqeq/tauchange'
%         Passed: 1
%         Failed: 0
%     Incomplete: 0
%       Duration: 1.7105e-02
% Totals:
%    1 Passed, 0 Failed, 0 Incomplete.
%    0.017105 seconds testing time.