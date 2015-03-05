%UT_MEANMCBERNOULLI_G  unit tests for meanMCBer_g and Test_meanMCBer_g 
classdef ut_meanMCBer_g < matlab.unittest.TestCase
    
    methods(Test)
        
        function meanMCBer_gOfabs(testCase)
            in_param.abstol=1e-3;
            in_param.alpha = 0.01;
            p=1/90;
            Yrand=@(n) rand(n,1)<p;
            pHat = meanMCBer_g(Yrand,in_param);
            actualerr = abs(pHat-p);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function meanMCBer_gOfrel(testCase)
            p=1/90;
            Yrand=@(n) rand(n,1)<p;
            in_param.abstol = 5e-2;
            pHat = meanMCBer_g(Yrand,'abstol',in_param.abstol);
            actualerr = abs(pHat-p);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function meanMCBer_gOfparsing(testCase)
            p=1/90;
            in_param.abstol = -1e-2;
            pHat = testCase.verifyWarning(@()meanMCBer_g...
                (@(n) (rand(n,1)<p).^2,...
                in_param),'MATLAB:meanMCBer_g:abstolneg');
            actualerr = abs(pHat-p);
            testCase.verifyLessThanOrEqual(actualerr,abs(in_param.abstol));
        end
        
        function meanMCBer_gOfWorkouts(testCase)
            [ut_abserr,ut_abstol] = Test_meanMCBer_g;
            testCase.verifyLessThanOrEqual(ut_abserr,ut_abstol);
        end
    end
end
