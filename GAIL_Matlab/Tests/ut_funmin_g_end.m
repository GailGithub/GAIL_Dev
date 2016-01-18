%UT_FUNMIN_G_end unit test for funmin_g at end points
classdef ut_funmin_g_end < matlab.unittest.TestCase
    
    methods(Test)
           
        function funmin_gEXM1(testCase)
            % minimum appears at left-end point
            f=@(x) sin(x);
            [fmin,out]=funmin_g(f);
            actualminima=0;
            testCase.verifyEqual(actualminima,out.intervals(1,1));
        end
        
        function funmin_gEXM2(testCase)
            % minimum appears at left-end point
            f=@(x) sin(1-x);
            [fmin,out]=funmin_g(f);
            actualminima=1;
            testCase.verifyEqual(actualminima,out.intervals(2,end));
        end
        
        function funmin_gEXM3(testCase)
            % minimum appears at both left and right-end points
            f=@(x) x.*(1-x);
            [fmin,out]=funmin_g(f);
            actualminima1=0;
            actualminima2=1;
            testCase.verifyEqual(actualminima1,out.intervals(1,1));
            testCase.verifyEqual(actualminima2,out.intervals(2,end));
        end
        
        function funmin_gEXM4(testCase)
            % minimum appears at both left- and right-end points
            f=@(x) sin(pi*x);
            [fmin,out]=funmin_g(f);
            actualminima1=0;
            actualminima2=1;
            testCase.verifyEqual(actualminima1,out.intervals(1,1));
            testCase.verifyEqual(actualminima2,out.intervals(2,end));
        end
        
        function funmin_gEXM5(testCase)
            % minimum appears at left end point
            f=@(x) -(x-1).^2+1;
            tol = 1e-7;
            Xtol = 1e-4;
            [fmin,out] = funmin_g(f,0,1,tol,Xtol,10,10,1000000);
            xmin_true = 0; 
            fmin_true = f(xmin_true);
            ferror = abs(fmin - fmin_true);
            xerror = abs(mean(out.intervals) - xmin_true);
            check = ferror < tol || xerror < Xtol;
            testCase.verifyEqual(check, true);
        end
        
        function funmin_gEXM6(testCase)
            % minimum appears at left end point
            f=@(x) -(x-0.3).^2+1; 
            tol = 1e-7;
            Xtol = 1e-4; 
            [fmin,out_param] = funmin_g(f,-2,2,tol,Xtol,10,10,1000000);
            xmin_true = -2; 
            fmin_true = f(xmin_true);
            ferror = abs(fmin - fmin_true);
            xerror = abs(mean(out_param.intervals) - xmin_true);
            check = ferror < tol || xerror < Xtol;
            testCase.verifyEqual(check, true);
        end
    end
end

