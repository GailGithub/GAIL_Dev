%UT_TRAUB_PAPER unit tests for the Traub paper
classdef ut_traub_paper < matlab.unittest.TestCase

    methods(Test)
        
        function test_traub_paper_funmin_g(testCase)
            nrep = 100; abstol = 1e-6; ninit = 250;
   
            [timeratio,npointsratio,matfilename1] = ...
                funappx_g_test(nrep,abstol,ninit,'funappxPenalty_g');
            
            [timeratio,npointsratio,matfilename2] = ...
                funmin_g_test(nrep,abstol,'funmin_g');

            [funappx_g_success_rate,funmin_g_success_rate] = ...
                LocallyAdaptivePaperFigs(matfilename1,matfilename2);
            
            testCase.verifyGreaterThanOrEqual(funappx_g_success_rate, [99 99 99]);
            testCase.verifyGreaterThanOrEqual(funmin_g_success_rate, [99 99 99]);
        end
    
  end
end