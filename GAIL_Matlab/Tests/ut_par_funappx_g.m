%UT_PAR_FUNAPPX_G unit tests for par_funappx_g
classdef ut_par_funappx_g < matlab.unittest.TestCase
  
  methods(Test)
    function par_funappx_gofConstantFunction(testCase)
      f = @(x) 3;
      in_param.maxiter = 1;
      in_param.nlo = 3;
      in_param.nhi = 3;
      workers = 8;
      [fappx, result] = par_funappx_g(workers,f,in_param);
      testCase.verifyLessThanOrEqual(result.iter, 1);
      testCase.verifyLessThanOrEqual(result.npoints, 5*workers-workers+1);
      x = 0:0.1:1;
      testCase.verifyLessThanOrEqual(norm(fappx(x)-f(x)), eps);
    end
    
    function par_funappx_gOfx(testCase)
      f = @(x) x;
      in_param.abstol = 10^(-8);
      in_param.nlo = 10;
      in_param.nhi = 100;
      workers = 8;
      [fappx, result] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.iter, 1);
    end
    
    function par_funappx_gOf100000x(testCase)
      f = @(x) 100000 .* x;
      in_param.abstol = 10^(-8);
      in_param.nlo = 5;
      in_param.nhi = 5;
      workers = 8;
      [fappx, result] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.iter, 1);
      testCase.verifyLessThanOrEqual(result.npoints, 5*workers-workers+1);
    end
    
    function par_funappx_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 10^(-8);
      in_param.nlo = 10;
      in_param.nhi = 100;
      workers = 8;
      [fappx, ~] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8);
      in_param.nlo = 10;
      in_param.nhi = 10;
      workers = 8;
      [fappx, ~] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfexponential(testCase)
      f = @(x) exp(-100*(x-0.7).^2);
      in_param.abstol = 10^(-8);
      in_param.nlo = 10;
      in_param.nhi = 100;
      workers = 8;
      [fappx, ~] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1);
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfxab(testCase)
      f = @(x) x;
      in_param.a = 0;
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      workers = 8;
      [fappx,result] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfxsquareab(testCase)
      f = @(x) x.^2;
      in_param.a = 0;
      in_param.b = 0.001;
      in_param.abstol = 10^(-8);
      workers = 8;
      [fappx, result] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfsinab(testCase)
      f = @(x) sin(x);
      in_param.a = 0;
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.nlo = 10;
      in_param.nhi = 100;
      workers = 8;
      [fappx, result] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfexponentialab(testCase)
      f = @(x)  exp(-100*(x-0.7).^2);
      in_param.a = 0;
      in_param.b = 10;
      in_param.abstol = 10^(-6);
      in_param.nlo = 10;
      in_param.nhi = 100;
      workers = 8;
      [fappx, result] = par_funappx_g(workers,f,in_param);
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function par_funappx_gOfblea(testCase)
      f = @(x) x.^2;
      in_param.a = 2;
      in_param.b = 1;
      workers = 8;
      [fappx, result] = testCase.verifyWarning(@()par_funappx_g(workers,f,in_param),'GAIL:gail1D_in_param:blea');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function par_funappx_gOfnofunction(testCase)
      workers = 8;
      [fappx, result] = testCase.verifyWarning(@()par_funappx_g(workers),'GAIL:gail_in_param:notfunction');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-result.f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function par_funappx_gOfbeqa(testCase)
      f = @(x) x.^2;
      in_param.a = 1;
      in_param.b = 1;
      workers = 8;
      [fappx, result] = testCase.verifyWarning(@()par_funappx_g(workers,f,in_param),'GAIL:gail1D_in_param:beqa');
      x = rand(10000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
    function par_funappx_gOfexceedbudget(testCase)
      f = @(x) x.^2;
      in_param.nmax = 200;
      workers = 8;
      %lastwarn('')
      [~, result] = ... %testCase.verifyWarning(@()
      par_funappx_g(workers,f,in_param);
      %,'GAIL:funappx_g:exceedbudget');
      %testCase.verifyEqual(isempty(lastwarn), false)   
      testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
      testCase.verifyEqual(result.exitflag,logical([1 0 0 0 0]));
    end
    
    function par_funappx_gOfexceediter(testCase)
      f = @(x) x.^2;
      in_param.maxiter = 2;
      workers = 8;
      %lastwarn('')
      [~, result] = ... %testCase.verifyWarning(@()
        par_funappx_g(workers,f,in_param);
      %,'GAIL:par_funappx_g:exceediter');
      %testCase.verifyEqual(isempty(lastwarn), 0);
      testCase.verifyEqual(result.maxiter,result.iter);
    end
    
    function par_funappx_gOnpointsoflinear(testCase)
      f = @(x) 3*x + 5;
      in_param.nlo = 5;
      in_param.nhi = 5;
      workers = 8;
      [~, result] = par_funappx_g(workers,f,in_param);
      testCase.verifyEqual(result.npoints,5*workers-workers+1);
    end
    
    function par_funappx_gOnpointsofconstant(testCase)
      f = @(x) 5;
      in_param.nlo = 1;
      in_param.nhi = 1;
      workers = 8;
      [~, result] = par_funappx_g(workers,f,in_param);
      testCase.verifyEqual(result.npoints,5*workers-workers+1);
    end
    
    function par_funappx_gHighAccuracy(testCase)
      workers = 8;
      [fappx, result] = par_funappx_g(workers,@(x)exp(x),'a',-2,'b',2,'nhi',20,'nlo',10, 'abstol', 1e-12);
      x = rand(1000000,1)*(result.b-result.a)+result.a;
      actualerr = max(abs(fappx(x)-result.f(x)));
      testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
  end
end
