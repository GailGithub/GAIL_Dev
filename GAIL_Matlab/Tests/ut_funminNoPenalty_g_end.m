%UT_FUNMINNOPENALTY_G_end unit test for funminNoPenalty_g at end points
classdef ut_funminNoPenalty_g_end < matlab.unittest.TestCase
    
    methods(Test)
           
        function funminNoPenalty_gEXM1(testCase)
            % minimum appears at left-end point
            f=@(x) sin(x);
            [fmin,out]=funminNoPenalty_g(f);
            actualminima=0;
            testCase.verifyEqual(actualminima,out.intervals(1,1));
        end
        
        function funminNoPenalty_gEXM2(testCase)
            % minimum appears at left-end point
            f=@(x) sin(1-x);
            [fmin,out]=funminNoPenalty_g(f);
            actualminima=1;
            testCase.verifyEqual(actualminima,out.intervals(2,end));
        end
        
        function funminNoPenalty_gEXM3(testCase)
            % minimum appears at both left and right-end points
            f=@(x) x.*(1-x);
            [fmin,out]=funminNoPenalty_g(f);
            actualminima1=0;
            actualminima2=1;
            testCase.verifyEqual(actualminima1,out.intervals(1,1));
            testCase.verifyEqual(actualminima2,out.intervals(2,end));
        end
        
        function funminNoPenalty_gEXM4(testCase)
            % minimum appears at both left- and right-end points
            f=@(x) sin(pi*x);
            [fmin,out]=funminNoPenalty_g(f);
            actualminima1=0;
            actualminima2=1;
            testCase.verifyEqual(actualminima1,out.intervals(1,1));
            testCase.verifyEqual(actualminima2,out.intervals(2,end));
        end
        
        function funminNoPenalty_gEXM5(testCase)
            % minimum appears at left end point
            f=@(x) -(x-1).^2+1;
            tol = 1e-7;
            Xtol = 1e-4;
            [fmin,out] = funminNoPenalty_g(f,0,1,tol,Xtol,10,10,1000000);
            xmin_true = 0; 
            fmin_true = f(xmin_true);
            ferror = abs(fmin - fmin_true);
            xerror = abs(mean(out.intervals) - xmin_true);
            check = ferror < tol || xerror < Xtol;
            testCase.verifyEqual(check, true);
        end
        
        function funminNoPenalty_gEXM6(testCase)
            % minimum appears at left end point
            f=@(x) -(x-0.3).^2+1; 
            tol = 1e-7;
            Xtol = 1e-4; 
            [fmin,out] = funminNoPenalty_g(f,-2,2,tol,Xtol,10,10,1000000);
            xmin_true = -2; 
            fmin_true = f(xmin_true);
            ferror = abs(fmin - fmin_true);
            xerror = abs(mean(out.intervals) - xmin_true);
            check = ferror < tol || xerror < Xtol;
            testCase.verifyEqual(check, true);
        end
    end
end

