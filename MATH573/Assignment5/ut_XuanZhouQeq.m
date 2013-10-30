classdef ut_XuanZhouQeq < matlab.unittest.TestCase
    methods(Test)
       function ut(testCase)
            testCase.verifyError(@()XuanZhouQeq('1','2','3'),'XuanZhouQeq:InputMustBeNumeric');
        end
    end
end

% run(ut_XuanZhouQeq)
% Running ut_XuanZhouQeq
% .
% Done ut_XuanZhouQeq
% __________
% 
% 
% ans = 
% 
%   TestResult with properties:
% 
%           Name: 'ut_XuanZhouQeq/ut'
%         Passed: 1
%         Failed: 0
%     Incomplete: 0
%       Duration: 0.0074
% 
% Totals:
%    1 Passed, 0 Failed, 0 Incomplete.
%    0.007401 seconds testing time.