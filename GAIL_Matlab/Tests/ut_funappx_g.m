%UT_FUNAPPX_G fast unit tests for funappx_g
classdef ut_funappx_g < matlab.unittest.TestCase

  methods(Test)
             
    function funappx_gOfx(testCase)
      f = @(x) x;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 100;
      [fappx, ~] = funappx_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 100;
      [fappx, ~] = funappx_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      [fappx, ~] = funappx_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfexponential(testCase)
      f = @(x) exp(-100*(x-0.7).^2);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 100;
      [fappx, ~] = funappx_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
          
    function funappx_gOfxab(testCase)
      f = @(x) x;
      in_param.a = 0; 
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      [fappx,result] = funappx_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfxsquareab(testCase)
      f = @(x) x.^2;
      in_param.a = 0; 
      in_param.b = 0.001;
      in_param.abstol = 10^(-8);
      [fappx, result] = funappx_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfsinab(testCase)
      f = @(x) sin(x);
      in_param.a = 0; 
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.nlo = 10;
      in_param.nhi = 100;
      [fappx, result] = funappx_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfexponentialab(testCase)
      f = @(x)  exp(-100*(x-0.7).^2);
      in_param.a = 0; 
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.nlo = 10;
      in_param.nhi = 100;
      [fappx, result] = funappx_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfblea(testCase)
      f = @(x) x.^2;
      in_param.a = 2; 
      in_param.b = 1;  
      [fappx, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:blea');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function funappx_gOfnofunction(testCase)
      [fappx, result] = testCase.verifyWarning(@()funappx_g,'MATLAB:funappx_g:nofunction');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-result.f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function funappx_gOfbeqa(testCase)
      f = @(x) x.^2;
      in_param.a = 1; 
      in_param.b = 1;  
      [fappx, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:beqa');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function funappx_gOfexceedbudget(testCase)
      f = @(x) x.^2;
      in_param.nmax = 1000; 
      [~, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:exceedbudget');
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
    end
    
    function funappx_gOfexceediter(testCase)
      f = @(x) x.^2;
      in_param.maxiter = 2;
      [~, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'MATLAB:funappx_g:exceediter');
      testCase.verifyLessThan(result.maxiter,result.iter);
    end
    
  end
end
