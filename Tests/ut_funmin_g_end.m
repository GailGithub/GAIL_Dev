%UT_FUNMIN_G_END unit test for funmin_g at end points
classdef ut_funmin_g_end < matlab.unittest.TestCase
    
    methods(Test)
        
        function funmin_gEXM1(testCase)
            % minimum appears at left end point
            f=@(x) -(x-1).^2+1;
            tol = 1e-7;
            [fmin,out] = funmin_g(f,0,1,tol,10,1000000);
            xmin_true = 0; 
            fmin_true = f(xmin_true);
            ferror = abs(fmin - fmin_true);
            check = ferror < tol;
            testCase.verifyEqual(check, true);
        end
        
        function funmin_gEXM2(testCase)
            % minimum appears at right end point
            f=@(x) -(x+0.3).^2+1; 
            tol = 1e-7;
            [fmin,out] = funmin_g(f,-2,2,tol,10,1000000);
            xmin_true = 2; 
            fmin_true = f(xmin_true);
            ferror = abs(fmin - fmin_true);
            check = ferror < tol;
            testCase.verifyEqual(check, true);
        end
    end
end


