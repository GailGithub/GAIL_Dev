%UT_MEANMCabs_G  unit test for meanMCabs_g
classdef ut_meanMCabs_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function meanMCabs_gOfexp(testCase)
      Y = @(n) exp(rand(n,1));
      in_param.abstol = 1e-2;
      meanY = meanMCabs_g(Y,in_param);
      exactY = exp(1)-1;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMCabs_gOfxsquare(testCase)
      Y = @(n) rand(n,1).^2;
      in_param.abstol = 1e-2;
      meanY = meanMCabs_g(Y,in_param);
      exactY = 1/3;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMCabs_gOfsin(testCase)
      Y = @(n) sin(rand(n,1));
      in_param.abstol = 1e-2;
      meanY = meanMCabs_g(Y,in_param);
      exactY = 1-cos(1);
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMCabs_gOfparsing(testCase)
      in_param.abstol = -1e-2;
      meanY = testCase.verifyWarning(@()meanMCabs_g(@(n) rand(n,1).^2,...
        in_param.abstol),'GAIL:meanMCabs_g:abstolneg');
      exactY = 1/3;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,abs(in_param.abstol));
    end
    
    function meanMC_gOfnonRandomInput(testCase)
        in_param.abstol = 1e-2;
        meanY = testCase.verifyWarning(@()meanMCabs_g(@(x) x.^2,...
            in_param.abstol),'GAIL:meanMCabs_g:yrandnotlengthN');
    end
  end
end

