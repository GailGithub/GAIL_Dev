%UT_funappxNoPenalty_g fast unit tests for funappxNoPenalty_g
classdef ut_funappxNoPenalty_g < matlab.unittest.TestCase
    
    methods(Test)
        function funappxNoPenalty_gofConstantFunction(testCase)
            f = @(x) 3;
            in_param.maxiter = 1;
            in_param.nlo = 3;
            in_param.nhi = 3;
            [fappx, result] = funappxNoPenalty_g(f,in_param);
            testCase.verifyLessThanOrEqual(result.iter, 1);
            testCase.verifyEqual(result.npoints,3);
            x = 0:0.1:1;
            testCase.verifyLessThanOrEqual(norm(fappx(x)-f(x)), eps);
        end
        
        function funappxNoPenalty_gOfx(testCase)
            f = @(x) x;
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, result] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.iter, 1); 
        end 
        
        function funappxNoPenalty_gOf100000x(testCase)
            f = @(x) 100000 .* x;
            in_param.abstol = 10^(-8);
            in_param.nlo = 5;
            in_param.nhi = 5;
            [fappx, result] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.iter, 1);
            testCase.verifyLessThanOrEqual(result.npoints, 5);
         end
        
        function funappxNoPenalty_gOfxsquare(testCase)
            f = @(x) x.^2;
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, ~] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfsin(testCase)
            f = @(x) sin(x);
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 10;
            [fappx, ~] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfexponential(testCase)
            f = @(x) exp(-100*(x-0.7).^2);
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, ~] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfxab(testCase)
            f = @(x) x;
            in_param.a = 0;
            in_param.b = 10;
            in_param.abstol = 10^(-6);
            [fappx,result] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfxsquareab(testCase)
            f = @(x) x.^2;
            in_param.a = 0;
            in_param.b = 0.001;
            in_param.abstol = 10^(-8);
            [fappx, result] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfsinab(testCase)
            f = @(x) sin(x);
            in_param.a = 0;
            in_param.b = 10;
            in_param.abstol = 10^(-6);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, result] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfexponentialab(testCase)
            f = @(x)  exp(-100*(x-0.7).^2);
            in_param.a = 0;
            in_param.b = 10;
            in_param.abstol = 10^(-6);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, result] = funappxNoPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxNoPenalty_gOfblea(testCase)
            f = @(x) x.^2;
            in_param.a = 2;
            in_param.b = 1;
            [fappx, result] = testCase.verifyWarning(@()funappxNoPenalty_g(f,in_param),'GAIL:funappxNoPenalty_g:blea');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappxNoPenalty_gOfnofunction(testCase)
            [fappx, result] = testCase.verifyWarning(@()funappxNoPenalty_g,'GAIL:funappxNoPenalty_g:nofunction');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-result.f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappxNoPenalty_gOfbeqa(testCase)
            f = @(x) x.^2;
            in_param.a = 1;
            in_param.b = 1;
            [fappx, result] = testCase.verifyWarning(@()funappxNoPenalty_g(f,in_param),'GAIL:funappxNoPenalty_g:beqa');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappxNoPenalty_gOfexceedbudget(testCase)
            f = @(x) x.^2;
            in_param.nmax = 200;
            [~, result] = testCase.verifyWarning(@()funappxNoPenalty_g(f,in_param),'GAIL:funappxNoPenalty_g:exceedbudget');
            testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
            testCase.verifyEqual(result.exit,([1 0]==1));
        end
        
        function funappxNoPenalty_gOfexceediter(testCase)
            f = @(x) x.^2;
            in_param.maxiter = 2;
            [~, result] = testCase.verifyWarning(@()funappxNoPenalty_g(f,in_param),'GAIL:funappxNoPenalty_g:exceediter');
            testCase.verifyEqual(result.maxiter,result.iter);
        end
        
         function funappxNoPenalty_gOnpointsoflinear(testCase)
            f = @(x) 3*x + 5;
            in_param.nlo = 5; in_param.nhi =5;
            [~, result] = funappxNoPenalty_g(f,in_param);
            testCase.verifyEqual(result.npoints,5);
         end
        
         function funappxNoPenalty_gOnpointsofconstant(testCase)
            f = @(x) 5;
            in_param.nlo = 5; in_param.nhi =5;
            [~, result] = funappxNoPenalty_g(f,in_param);
            testCase.verifyEqual(result.npoints,5);
         end
         
    end
end
