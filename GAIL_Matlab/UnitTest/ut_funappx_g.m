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
      x = rand(1000,1);
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1);
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1);
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1);
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-f(x)));
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
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappx_gOfblea(testCase)
      f = @(x) x.^2;
      in_param.a = 2; 
      in_param.b = 1;  
      [appxf, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:blea');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
    function funappx_gOfnofunction(testCase)
      [appxf, result] = testCase.verifyWarning(@()funappx_g,'MATLAB:funappx_g:nofunction');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-result.f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
    function funappx_gOfbeqa(testCase)
      f = @(x) x.^2;
      in_param.a = 1; 
      in_param.b = 1;  
      [appxf, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:beqa');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(appxf(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
  end
end
