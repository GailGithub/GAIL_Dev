%UT_funappx_g fast unit tests for funappx_g
classdef ut_funappx_g < matlab.unittest.TestCase
    
    methods(Test)
        function funappx_gofConstantFunction(testCase)
            f = @(x) 3;
            in_param.maxiter = 1;
            in_param.ninit = 5;
            [fappx, result] = funappx_g(f,in_param);
            testCase.verifyLessThanOrEqual(result.iter, 1);
            testCase.verifyEqual(result.npoints, 6);
            x = 0:0.1:1;
            testCase.verifyLessThanOrEqual(norm(fappx(x)-f(x)), eps);
        end
        
        function funappx_gOfx(testCase)
            f = @(x) x;
            in_param.abstol = 10^(-8);
            in_param.ninit = 35;
            [fappx, result] = funappx_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.iter, 1); 
        end 
        
        function funappx_gOf100000x(testCase)
            f = @(x) 100000 .* x;
            in_param.abstol = 10^(-8);
            in_param.ninit = 5;
            [fappx, result] = funappx_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.iter, 1);
            testCase.verifyLessThanOrEqual(result.npoints, 6);
         end
        
        function funappx_gOfxsquare(testCase)
            f = @(x) x.^2;
            in_param.abstol = 10^(-8);
            in_param.ninit = 35;
            [fappx, ~] = funappx_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappx_gOfsin(testCase)
            f = @(x) sin(x);
            in_param.abstol = 10^(-8);
            in_param.ninit = 10;
            [fappx, ~] = funappx_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappx_gOfexponential(testCase)
            f = @(x) exp(-100*(x-0.7).^2);
            in_param.abstol = 10^(-8);
            in_param.ninit = 35;
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
            in_param.ninit = 35;
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
            [fappx, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'GAIL:funappx_g:blea');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappx_gOfnofunction(testCase)
            [fappx, result] = testCase.verifyWarning(@()funappx_g,'GAIL:funappx_g:nofunction');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-result.f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappx_gOfbeqa(testCase)
            f = @(x) x.^2;
            in_param.a = 1;
            in_param.b = 1;
            [fappx, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'GAIL:funappx_g:beqa');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappx_gOfexceedbudget(testCase)
            f = @(x) x.^2;
            in_param.nmax = 200;
            [~, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'GAIL:funappx_g:exceedbudget');
            testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
            testCase.verifyEqual(result.exitflag,logical([1 0]));
        end
        
        function funappx_gOfexceediter(testCase)
            f = @(x) x.^2;
            in_param.maxiter = 2;
            [~, result] = testCase.verifyWarning(@()funappx_g(f,in_param),'GAIL:funappx_g:exceediter');
            testCase.verifyEqual(result.maxiter,result.iter);
        end
        
         function funappx_gOnpointsoflinear(testCase)
            f = @(x) 3*x + 5;
            in_param.ninit = 5;
            [~, result] = funappx_g(f,in_param);
            testCase.verifyEqual(result.npoints,6);
         end
        
         function funappx_gOnpointsofconstant(testCase)
            f = @(x) 5;
            in_param.ninit = 5;
            [~, result] = funappx_g(f,in_param);
            testCase.verifyEqual(result.npoints,6);
         end
         
         
    end
end

