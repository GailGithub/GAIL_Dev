%UT_FUNAPPX_G unit test for funappxab_g
classdef ut_funappx_g < matlab.unittest.TestCase

  methods(Test)
    
         
    function funappx_gOfx(testCase)
      f = @(x) x;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      in_param.nmax = 10^6;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      in_param.nmax = 10^6;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      in_param.nmax = 10^6;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfexponential(testCase)
      f = @(x) exp(-100*(x-0.7).^2);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 100;
      in_param.nhi = 100;
      in_param.nmax = 10^6;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)-1;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
          
    function funappx_gOfxab(testCase)
      f = @(x) x;
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^6;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfxsquareab(testCase)
      f = @(x) x.^2;
      in_param.a = 0; 
      in_param.b = 0.001;
      in_param.abstol = 10^(-8);
      in_param.nmax = 10^8;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfsinab(testCase)
      f = @(x) sin(x);
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^8;
      in_param.nlo = 100;
      in_param.nhi = 1000;
      [appxf, result] = funappx_g(f,in_param);
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfexponentialab(testCase)
      f = @(x)  exp(-100*(x-0.7).^2);
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^8;
      in_param.nlo = 100;
      in_param.nhi = 1000;
      [appxf, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:peaky');
      x = sqrt(2)*in_param.b/2;
      actualerr = abs(appxf(x)-f(x));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
  end
end
