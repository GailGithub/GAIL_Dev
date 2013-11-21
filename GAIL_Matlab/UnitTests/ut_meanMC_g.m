%UT_MEANMC_G unit test for meanMC_g
classdef ut_meanMC_g < matlab.unittest.TestCase
    
    methods(Test)
        
        function meanMC_gOfexp(testCase)
            Y = @(n) exp(rand(n,1));
            in_param.abstol = 1e-2;
            meanY = meanMC_g(Y,in_param);
            exactY = exp(1)-1;
            actualerr = abs(meanY-exactY);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function meanMC_gOfxsquare(testCase)
            Y = @(n) rand(n,1).^2;
            in_param.abstol = 1e-2;
            meanY = meanMC_g(Y,in_param);
            exactY = 1/3;
            actualerr = abs(meanY-exactY);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function meanMC_gOfsin(testCase)
            Y = @(n) sin(rand(n,1));
            in_param.abstol = 1e-2;
            meanY = meanMC_g(Y,in_param);
            exactY = 1-cos(1);
            actualerr = abs(meanY-exactY);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function meanMC_gOfparsing(testCase)
            in_param.abstol = -1e-2;
            meanY = testCase.verifyWarning(@()meanMC_g(@(n) rand(n,1).^2,...
                in_param.abstol),'MATLAB:meanMC_g:abstolneg');
            exactY = 1/3;
            actualerr = abs(meanY-exactY);
            testCase.verifyLessThanOrEqual(actualerr,abs(in_param.abstol));
        end
    end
end
