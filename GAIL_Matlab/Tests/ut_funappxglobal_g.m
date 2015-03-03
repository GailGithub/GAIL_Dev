%UT_FUNAPPXGLOBAL_G fast unit tests for funappxglobal_g
classdef ut_funappxglobal_g < matlab.unittest.TestCase

  methods(Test)
             
    function funappxglobal_gOfx(testCase)
      f = @(x) x;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      in_param.nmax = 10^6;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      in_param.nmax = 10^6;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      in_param.nmax = 10^6;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfexponential(testCase)
      f = @(x) exp(-100*(x-0.7).^2);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 100;
      in_param.nhi = 100;
      in_param.nmax = 10^6;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
          
    function funappxglobal_gOfxab(testCase)
      f = @(x) x;
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^6;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfxsquareab(testCase)
      f = @(x) x.^2;
      in_param.a = 0; 
      in_param.b = 0.001;
      in_param.abstol = 10^(-8);
      in_param.nmax = 10^8;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfsinab(testCase)
      f = @(x) sin(x);
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^8;
      in_param.nlo = 100;
      in_param.nhi = 1000;
      [pp, result] = funappxglobal_g(f,in_param);
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfexponentialab(testCase)
      f = @(x)  exp(-100*(x-0.7).^2);
      in_param.a = 0; 
      in_param.b = 1000;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^8;
      in_param.nlo = 100;
      in_param.nhi = 1000;
      [pp, result] = testCase.verifyWarning(@()funappxglobal_g(f,in_param),'MATLAB:funappxglobal_g:peaky');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funappxglobal_gOfblea(testCase)
      f = @(x) x.^2;
      in_param.a = 2; 
      in_param.b = 1;  
      [pp, result] = testCase.verifyWarning(@()funappxglobal_g(f,in_param),'MATLAB:funappxglobal_g:blea');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
    function funappxglobal_gOfnofunction(testCase)
      [pp, result] = testCase.verifyWarning(@()funappxglobal_g,'MATLAB:funappxglobal_g:nofunction');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-result.f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
    function funappxglobal_gOfbeqa(testCase)
      f = @(x) x.^2;
      in_param.a = 1; 
      in_param.b = 1;  
      [pp, result] = testCase.verifyWarning(@()funappxglobal_g(f,in_param),'MATLAB:funappxglobal_g:beqa');
      x = rand(1000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
  end
end
