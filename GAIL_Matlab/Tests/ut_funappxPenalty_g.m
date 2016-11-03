%UT_FUNAPPX_G fast unit tests for funappxPenalty_g
classdef ut_funappxPenalty_g < matlab.unittest.TestCase
    
    methods(Test)
        function funappxPenalty_gofConstantFunction(testCase)
            f = @(x) 3;
            in_param.maxiter = 1;
            in_param.nlo = 1;
            in_param.nhi = 1;
            [fappx, result] = testCase.verifyWarning(@()funappxPenalty_g(f,in_param),'GAIL:funappxPenalty_g:exceediter');
            testCase.verifyLessThanOrEqual(result.iter, 1);
            testCase.verifyEqual(result.npoints,3);
            x = 0:0.1:1;
            testCase.verifyLessThanOrEqual(norm(fappx(x)-f(x)), eps);
        end
        
        function funappxPenalty_gOfx(testCase)
            f = @(x) x;
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, result] = funappxPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.iter, 1); 
        end 
        
        function funappxPenalty_gOf100000x(testCase)
            f = @(x) 100000 .* x;
            in_param.abstol = 10^(-8);
            in_param.nlo = 1;
            in_param.nhi = 1;
            [fappx, result] = funappxPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.iter, 1);
            testCase.verifyLessThanOrEqual(result.npoints, 3);
         end
        
        function funappxPenalty_gOfxsquare(testCase)
            f = @(x) x.^2;
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, ~] = funappxPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfsin(testCase)
            f = @(x) sin(x);
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 10;
            [fappx, ~] = funappxPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfexponential(testCase)
            f = @(x) exp(-100*(x-0.7).^2);
            in_param.abstol = 10^(-8);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, ~] = funappxPenalty_g(f,in_param);
            x = rand(10000,1);
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfxab(testCase)
            f = @(x) x;
            in_param.a = 0;
            in_param.b = 10;
            in_param.abstol = 10^(-6);
            [fappx,result] = funappxPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfxsquareab(testCase)
            f = @(x) x.^2;
            in_param.a = 0;
            in_param.b = 0.001;
            in_param.abstol = 10^(-8);
            [fappx, result] = funappxPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfsinab(testCase)
            f = @(x) sin(x);
            in_param.a = 0;
            in_param.b = 10;
            in_param.abstol = 10^(-6);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, result] = funappxPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfexponentialab(testCase)
            f = @(x)  exp(-100*(x-0.7).^2);
            in_param.a = 0;
            in_param.b = 10;
            in_param.abstol = 10^(-6);
            in_param.nlo = 10;
            in_param.nhi = 100;
            [fappx, result] = funappxPenalty_g(f,in_param);
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function funappxPenalty_gOfblea(testCase)
            f = @(x) x.^2;
            in_param.a = 2;
            in_param.b = 1;
            [fappx, result] = testCase.verifyWarning(@()funappxPenalty_g(f,in_param),'GAIL:funappxPenalty_g:blea');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappxPenalty_gOfnofunction(testCase)
            [fappx, result] = testCase.verifyWarning(@()funappxPenalty_g,'GAIL:funappxPenalty_g:nofunction');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-result.f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappxPenalty_gOfbeqa(testCase)
            f = @(x) x.^2;
            in_param.a = 1;
            in_param.b = 1;
            [fappx, result] = testCase.verifyWarning(@()funappxPenalty_g(f,in_param),'GAIL:funappxPenalty_g:beqa');
            x = rand(10000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
        
        function funappxPenalty_gOfexceedbudget(testCase)
            f = @(x) x.^2;
            in_param.nmax = 1000;
            [~, result] = testCase.verifyWarning(@()funappxPenalty_g(f,in_param),'GAIL:funappxPenalty_g:exceedbudget');
            testCase.verifyLessThanOrEqual(result.npoints,result.nmax);
        end
        
        function funappxPenalty_gOfexceediter(testCase)
            f = @(x) x.^2;
            in_param.maxiter = 2;
            [~, result] = testCase.verifyWarning(@()funappxPenalty_g(f,in_param),'GAIL:funappxPenalty_g:exceediter');
            testCase.verifyEqual(result.maxiter,result.iter);
        end
        
         function funappxPenalty_gOnpointsoflinear(testCase)
            f = @(x) 3*x + 5;
            in_param.nlo = 1; in_param.nhi =1;
            [~, result] = funappxPenalty_g(f,in_param);
            testCase.verifyEqual(result.npoints,3);
         end
        
         function funappxPenalty_gOnpointsofconstant(testCase)
            f = @(x) 5;
            in_param.nlo = 1; in_param.nhi =1;
            [~, result] = funappxPenalty_g(f,in_param);
            testCase.verifyEqual(result.npoints,3);
         end

                 
         function funappxPenalty_gHighAccuracy(testCase)
            [fappx, result] = funappxPenalty_g(@(x)exp(x),'a',-2,'b',2,'nhi',20,'nlo',10, 'abstol', 1e-12);
            x = rand(1000000,1)*(result.b-result.a)+result.a;
            actualerr = max(abs(fappx(x)-result.f(x)));
            testCase.verifyLessThanOrEqual(actualerr,result.abstol);
        end
    end
end
