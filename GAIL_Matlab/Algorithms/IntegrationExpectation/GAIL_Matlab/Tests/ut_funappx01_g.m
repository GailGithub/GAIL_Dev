%UT_FUNAPPX01_G fast unit tests for funappx01_g
classdef ut_funappx01_g < matlab.unittest.TestCase

  methods(Test)
    
      
    function funappx01_gOfx(testCase)
      f = @(x) x;
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [appxf, result] = funappx01_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx01_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [appxf, result] = funappx01_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx01_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [appxf, result] = funappx01_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx01_gOfexponential(testCase)
      f = @(x) exp(-100*(x-0.7).^2);
      in_param.abstol = 10^(-8); 
      in_param.ninit = 100;
      in_param.nmax = 10^6;
      [appxf, result] = funappx01_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
  end
end
