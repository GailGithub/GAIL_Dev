classdef ut_Xinqeq < matlab.unittest.TestCase
% unit test for Xinqeq
methods (Test)
    function testXinSolution1(testCase)
        abstol= 1e-5;
        XinSolution = Xinqeq(1,-3,2);
        expSolution = [1,2];
        error = max(abs(XinSolution-expSolution));
        testCase.verifyLessThanOrEqual(error,abstol);
    end
    
    function testXinSolution2(testCase)
        abstol= 1e-5;
        XinSolution = Xinqeq([0;-1;1]);
        expSolution = 1;
        error = max(abs(XinSolution-expSolution));
        testCase.verifyLessThanOrEqual(error,abstol);
    end

    function testXinSolution3(testCase)
        abstol= 1e-5;
        coef.a=1; coef.b=3; coef.c=2;
        XinSolution = Xinqeq(coef);
        expSolution = [-2,-1];
        error = max(abs(XinSolution-expSolution));
        testCase.verifyLessThanOrEqual(error,abstol);
    end
end
end

% Running ut_Xinqeq
% .Warning: The polynomial is linear. 
% > In Xinqeq>Xinqeq_param at 126
%   In Xinqeq at 55
%   In ut_Xinqeq>ut_Xinqeq.testXinSolution2 at 14
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
% ..
% Done ut_Xinqeq
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
%    0.023348 seconds testing time.

