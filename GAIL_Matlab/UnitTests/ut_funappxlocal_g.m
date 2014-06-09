%UT_FUNAPPXLOCAL_G unit test for funappxlocal_g
classdef ut_funappxlocal_g < matlab.unittest.TestCase

  methods(Test)
             
    function funappx_gOfx(testCase)
      f = @(x) x;
      in_param.abstol = 10^(-8); 
      in_param.taulo = 10;
      in_param.tauhi = 10;
      [pp, ~] = funappxlocal_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 10^(-8); 
      in_param.nlo = 10;
      in_param.nhi = 10;
      [pp, ~] = funappxlocal_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8); 
      in_param.taulo = 10;
      in_param.tauhi = 10;
      [pp, ~] = funappx_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfexponential(testCase)
      f = @(x) exp(-100*(x-0.7).^2);
      in_param.abstol = 10^(-8); 
      in_param.nlo = 100;
      in_param.nhi = 100;
      [pp, ~] = funappxlocal_g(f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
          
    function funappx_gOfxab(testCase)
      f = @(x) x;
      in_param.a = 0; 
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.nmax = 10^6;
      [pp,result] = funappxlocal_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfxsquareab(testCase)
      f = @(x) x.^2;
      in_param.a = 0; 
      in_param.b = 0.001;
      in_param.abstol = 10^(-8);
      [pp, result] = funappxlocal_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfsinab(testCase)
      f = @(x) sin(x);
      in_param.a = 0; 
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.taulo = 100;
      in_param.tauhi = 1000;
      [pp, result] = funappxlocal_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfexponentialab(testCase)
      f = @(x)  exp(-100*(x-0.7).^2);
      in_param.a = 0; 
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.taulo = 100;
      in_param.tauhi = 1000;
      [pp, result] = funappxlocal_g(f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function funappx_gOfblea(testCase)
      f = @(x) x.^2;
      in_param.a = 2; 
      in_param.b = 1;  
      [pp, result] = testCase.verifyWarning(@()funappxlocal_g(f,in_param),'MATLAB:funappx_g:blea');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function funappx_gOfnofunction(testCase)
      [pp, result] = testCase.verifyWarning(@()funappxlocal_g,'MATLAB:funappx_g:nofunction');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-result.f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function funappx_gOfbeqa(testCase)
      f = @(x) x.^2;
      in_param.a = 1; 
      in_param.b = 1;  
      [pp, result] = testCase.verifyWarning(@()funappxlocal_g(f,in_param),'MATLAB:funappx_g:beqa');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(ppval(pp,x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
  end
end
