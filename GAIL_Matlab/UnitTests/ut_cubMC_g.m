%UT_CUBMC_G unit test for cubMC_g
classdef ut_cubMC_g < matlab.unittest.TestCase

  methods(Test)
      
    function cubMC_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 1e-2; 
      interval=[0;1];
      meanf = cubMC_g(f,interval,in_param);
      exactf = 0.33;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMC_gOfexp(testCase)
      f = @(x) exp(x);
      in_param.abstol = 1e-3; 
      interval=[0;1];
      meanf = cubMC_g(f,interval,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMC_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 1e-3; 
      interval=[0;1];
      meanf = cubMC_g(f,interval,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end 
    
    function meanMC_gOfmultierrfun(testCase)
        f = @(x) exp(-x(:,1).^2-x(:,2).^2);
        in_param.abstol = 1e-3;
        interval=[0 0;1 1];
        meanf = cubMC_g(f,interval,in_param);
        exactf = pi/4*erf(1)^2;
        actualerr = abs(meanf-exactf);
        testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
  end
end
