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
            yuhansol = yuhanqeq(0,1,2);
            exactsol= -2;        
            testCase.verifyEqual(yuhansol,exactsol);
        end
        
    end
end

% Running ut_yuhanqeq
% ..Warning:  coefficient of quadratic term is zero; quadratic equation becomes a linear
% equation. 
% > In yuhanqeq>qeq_param at 169
%   In yuhanqeq at 96
%   In ut_yuhanqeq>ut_yuhanqeq.yuhanqeqof012 at 23
%   In TestRunner>TestRunner.evaluateMethod at 511
%   In TestRunner>TestRunner.runMethodSettingCurrentException at 504
%   In TestRunner>TestRunner.runMethod at 601
%   In TestRunnerPlugin>TestRunnerPlugin.runMethod at 88
%   In TestRunnerPlugin>TestRunnerPlugin.runMethod at 88
%   In TestRunner>TestRunner.runMethodOnPlugins at 463
%   In TestRunner>TestRunner.runSingleMethodOnTestContentAndHandleException at 487
%   In TestRunner>TestRunner.runAllMethods at 497
%   In TestRunner>TestRunner.runMethodsOnTestContent at 474
%   In TestRunner>TestRunner.runTestMethod at 652
%   In TestRunnerPlugin>TestRunnerPlugin.runTestMethod at 84
%   In TestRunnerPlugin>TestRunnerPlugin.runTestMethod at 84
%   In TestRunner>TestRunner.runMethodOnPlugins at 463
%   In TestRunner>TestRunner.runTestSuite at 578
%   In TestRunnerPlugin>TestRunnerPlugin.runTestSuite at 52
%   In TestRunnerPlugin>TestRunnerPlugin.runTestSuite at 52
%   In FailureDiagnosticsPlugin>FailureDiagnosticsPlugin.runTestSuite at 77
%   In TestRunner>TestRunner.runMethodOnPlugins at 463
%   In TestRunner>TestRunner.run at 171
%   In TestSuite>TestSuite.run at 311
%   In RunnableTestContent>RunnableTestContent.run at 43 
% .
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
%    0.023164 seconds testing time.