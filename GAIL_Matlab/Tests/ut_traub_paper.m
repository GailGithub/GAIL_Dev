%UT_TRAUB_PAPER unit tests for the Traub paper
classdef ut_traub_paper < matlab.unittest.TestCase

    methods(Test)
        
        function test_traub_paper_funappx_g(testCase)
            nrep = 1000; abstol = 1e-6; ninit = 250;
   
            [timeratio,npointsratio,matfilename] = ...
                funappx_g_test(nrep,abstol,ninit,'funappxPenalty_g');

            [funappx_g_success_rate] =  LocallyAdaptivePaperFigs(matfilename);
            
            testCase.verifyGreaterThanOrEqual(funappx_g_success_rate, [99 99 99]);
            
        end
    
        function test_traub_paper_funmin_g(testCase)
            nrep = 1000; abstol = 1e-6;  
            
            [timeratio,npointsratio,matfilename] = ...
                funmin_g_test(nrep,abstol,'funmin_g');
            
            [~, funmin_g_success_rate] = LocallyAdaptivePaperFigs([],matfilename);
            testCase.verifyGreaterThanOrEqual(funmin_g_success_rate, [99 99 99]);
        end
    
  end
end
