%ut_cubSobol_g  unit test for cubSobol_g
classdef ut_cubSobol_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubSobol_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 1e-2;
      hyperbox = [0;1];
      meanf = cubSobol_g(f,hyperbox,in_param);
      exactf = 0.33;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubSobol_gOfexp(testCase)
      f = @(x) exp(x);
      in_param.abstol = 1e-3;
      hyperbox = [0;1];
      meanf = cubSobol_g(f,hyperbox,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubSobol_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 1e-3;
      hyperbox = [0;1];
      meanf = cubSobol_g(f,hyperbox,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubSobol_gOfmultierrfun(testCase)
      f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      in_param.abstol = 1e-3;
      hyperbox = [0 0;1 1];
      meanf = cubSobol_g(f,hyperbox,in_param);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubSobol_gOfwarning(testCase)
        testCase.verifyWarning(@()cubSobol_g,'MATLAB:cubSobol_g:fdnotgiven');
    end
    function cubSobol_gOdwarning(testCase)
        testCase.verifyWarning(@()cubSobol_g(@(x)x.^2,1.5),'MATLAB:cubSobol_g:hyperbox_error1');
    end
  end
end
