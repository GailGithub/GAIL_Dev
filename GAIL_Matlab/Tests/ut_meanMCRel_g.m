%UT_MEANMC_G  unit test for meanMC_g
classdef ut_meanMCRel_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function meanMCRel_gOfexp(testCase)
      in_param.abstol = 1e-2;
      meanY = testCase.verifyWarning(@()meanMCRel_g( @(n) exp(rand(n,1)),...
        in_param.abstol),'MATLAB:meanMCRel_g:maxreached');
      exactY = exp(1)-1;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMCRel_gOfxsquare(testCase)
      in_param.reltol = 5e-2;
      meanY = testCase.verifyWarning(@()meanMCRel_g(@(n) rand(n,1).^2,...
        in_param.reltol),'MATLAB:meanMCRel_g:maxreached');
      exactY = 1/3;
      actualerr = abs(meanY-exactY)/exactY;
      testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);
    end
    
    function meanMCRel_gOfsin(testCase)
      in_param.reltol = 1e-2;
      meanY = testCase.verifyWarning(@()meanMCRel_g(@(n) sin(rand(n,1)),...
        in_param.reltol),'MATLAB:meanMCRel_g:maxreached');
    exactY = 1-cos(1);
      actualerr = abs(meanY-exactY)/exactY;
      testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);
    end
    
    function meanMCRel_gOfparsing(testCase)
      in_param.abstol = -1e-2;
      meanY = testCase.verifyWarning(@()meanMCRel_g(@(n) rand(n,1).^2,...
        in_param.abstol),'MATLAB:meanMCRel_g:abstolneg');
      exactY = 1/3;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,abs(in_param.abstol));
    end
    
    function meanMCRel_gOfnonRandomInput(testCase)
        in_param.abstol = 1e-2;
        testCase.verifyWarning(@()meanMCRel_g(@(x) x.^2,...
            in_param.abstol),'MATLAB:meanMCRel_g:yrandnotlengthN');
        
    end
  end
end
