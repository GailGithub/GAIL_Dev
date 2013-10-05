%UT_MEANMC_G unit test for meanMC_g
classdef ut_meanMC_g < matlab.unittest.TestCase

  methods(Test)
      
    function meanMC_gOfexp(testCase)
      f = @(n) exp(rand(n,1));
      in_param.abstol = 1e-2; 
      meanf = meanMC_g(f,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMC_gOfxsquare(testCase)
      f = @(n) rand(n,1).^2;
      in_param.abstol = 1e-2; 
      meanf = meanMC_g(f,in_param);
      exactf = 1/3;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMC_gOfsin(testCase)
      f = @(n) sin(rand(n,1));
      in_param.abstol = 1e-2; 
      meanf = meanMC_g(f,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMC_gOfparsing(testCase)
      in_param.abstol = 1e-2; 
      meanf = meanMC_g(@(n) rand(n,1).^2,-1e-2);
      exactf = 1/3;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
  end
end
