%UT_FUNAPPX_G unit test for funappxab_g
classdef ut_funappxab_g < matlab.unittest.TestCase

  methods(Test)
    
      
    function funappxab_gOfx(testCase)
      f = @(x) x;
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^6;
      [appxf, result] = funappxab_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxab_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.a = 0; 
      in_param.b = 0.001;
      in_param.abstol = 10^(-8);
      in_param.nmax = 10^8;
      [appxf, result] = funappxab_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxab_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^8;
      [appxf, result] = funappxab_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxab_gOfexponential(testCase)
      f = @(x) sin(x);
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^8;
      [appxf, result] = funappxab_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
  end
end
