%UT_FUNMIN_G_CSC unit test for funmin_g_CSC_CSC
classdef ut_funmin_g_CSC < matlab.unittest.TestCase
    
    methods(Test)
           
        function funmin_g_CSCEXM1(testCase)
            in_param.abstol = 10^(-8); 
            in_param.TolX = 0;
            in_param.nmax = 10^7; 
            a = 0.027;
            z = 0.86;
            f = @(x) 0.5/a^2.*(-4*a.^2-(x-z).^2-(x-z-a).*abs(x-z-a)...
                +(x-z+a).*abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);
            [fmin, result] = funmin_g_CSC(f,in_param); 
            actualmin = -1;
            actualerr = abs(actualmin-fmin);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_g_CSCEXM2(testCase)
            in_param.abstol = 0; 
            in_param.TolX = 10^(-5);
            in_param.nmax = 10^6; 
            a = 0.012;
            z = 0.33;
            f = @(x) 0.5/a^2.*(-4*a.^2-(x-z).^2-(x-z-a).*abs(x-z-a)...
                +(x-z+a).*abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);
            [fmin, result] = funmin_g_CSC(f,in_param); 
            intnum = size(result.intervals,2);
            exactsolu = z;
            testflag = 1;
            for k=1:intnum
                if exactsolu <= result.intervals(2,k) && exactsolu >= ...
                        result.intervals(1,k)
                 testflag = 0;
                end
            end
            testCase.verifyLessThan(testflag,1);
            testCase.verifyLessThanOrEqual(result.volumeX,in_param.TolX);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_g_CSCEXM3(testCase)
            in_param.abstol = 10^(-8); 
            in_param.TolX = 10^(-6);
            in_param.nmax = 10^7; 
            a = 0.019;
            z = 0.77;
            f = @(x) 0.5/a^2.*(-4*a.^2-(x-z).^2-(x-z-a).*abs(x-z-a)...
                +(x-z+a).*abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);
            [fmin, result] = funmin_g_CSC(f,in_param); 
            actualmin = -1;
            actualerr = abs(actualmin-fmin);
            intnum = size(result.intervals,2);
            exactsolu = z;
            testflag = 1;
            if actualerr <= in_param.abstol
                testflag = 0;
            end
            for k=1:intnum
                if exactsolu <= result.intervals(2,k) && exactsolu >= ...
                        result.intervals(1,k)
                 testflag = 0;
                end
            end
            testCase.verifyLessThan(testflag,1);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_g_CSCEXM4(testCase)
            in_param.abstol = 0; 
            in_param.TolX = 10^(-5);
            in_param.nmax = 10^6; 
            a1=5; b1=10; c1=0.5-0.5*rand(1,1);
            a2=1; b2=10; c2=0.5+0.5*rand(1,1);
            f=@(x) -a1*exp(-(b1*(x-c1)).^2)-a2*exp(-(b2*(x-c2)).^2);
            exactsolu = fminbnd(f,0,(c1+c2)/2,optimset('TolX',in_param.TolX*10^(-2)));
            [fmin, result] = funmin_g_CSC(f,in_param); 
            intnum = size(result.intervals,2);
            testflag = 1;
            for k=1:intnum
                if exactsolu <= result.intervals(2,k) && exactsolu >= ...
                        result.intervals(1,k)
                 testflag = 0;
                end
            end
            testCase.verifyLessThan(testflag,1);
            testCase.verifyLessThanOrEqual(result.volumeX,in_param.TolX);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_g_CSCEXM5(testCase)
            in_param.abstol = 0; 
            in_param.TolX = 10^(-5);
            in_param.nmax = 10^6; 
            f=@(x) cos(4*pi*(x-0.2));
            exactsolu1 = 0.45;
            exactsolu2 = 0.95;
            [fmin, result] = funmin_g_CSC(f,in_param); 
            intnum = size(result.intervals,2);
            testflag1 = 1;
            testflag2 = 1;
            for k=1:intnum
                if exactsolu1 <= result.intervals(2,k) && exactsolu1 >= ...
                        result.intervals(1,k)
                 testflag1 = 0;
                end
                if exactsolu2 <= result.intervals(2,k) && exactsolu2 >= ...
                        result.intervals(1,k)
                 testflag2 = 0;
                end
            end
            testflag = testflag1 || testflag2;
            testCase.verifyLessThan(testflag,1);
            testCase.verifyLessThanOrEqual(result.volumeX,in_param.TolX);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_g_CSCEXM6(testCase)
            in_param.abstol = 0; 
            in_param.TolX = 10^(-3);
            in_param.nmax = 10^6; 
            in_param.a = -5;
            in_param.b = 5;
            f=@(x) sin(pi*x);
            exactsolu1 = -4.5;
            exactsolu2 = -2.5;
            exactsolu3 = -0.5;
            exactsolu4 = 1.5;
            exactsolu5 = 3.5;
            [fmin, result] = funmin_g_CSC(f,in_param); 
            intnum = size(result.intervals,2);
            testflag = ones(5,1);
            for k=1:intnum
                if exactsolu1 <= result.intervals(2,k) && exactsolu1 >= ...
                        result.intervals(1,k)
                 testflag(1) = 0;
                end
                if exactsolu2 <= result.intervals(2,k) && exactsolu2 >= ...
                        result.intervals(1,k)
                 testflag(2) = 0;
                end
                if exactsolu3 <= result.intervals(2,k) && exactsolu3 >= ...
                        result.intervals(1,k)
                 testflag(3) = 0;
                end
                if exactsolu4 <= result.intervals(2,k) && exactsolu4 >= ...
                        result.intervals(1,k)
                 testflag(4) = 0;
                end
                if exactsolu5 <= result.intervals(2,k) && exactsolu5 >= ...
                        result.intervals(1,k)
                 testflag(5) = 0;
                end
            end
            testCase.verifyEqual(sum(testflag),0);
            testCase.verifyLessThanOrEqual(result.volumeX,in_param.TolX);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end  
        
        function funmin_g_CSCEXM7(testCase)
            in_param.abstol = 10^(-3); 
            in_param.TolX = 0;
            in_param.nmax = 10^7;  
            z = 0.5;
            f = @(x) exp(-1./((x-z).^2)); % flat bottom
            [fmin,result] = funmin_g_CSC(f,in_param);
            testCase.verifyLessThanOrEqual(result.errest,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
            v=sqrt(-1/(log(in_param.abstol)));
            exactsolu1 =[.5-v, .5+v];
            % verify guarantee |min(f)-fmin| <= abstol
            testCase.verifyLessThanOrEqual(abs(fmin-0), in_param.abstol)
            % verify result.intervals is a subinterval of exactsolu1
            testCase.verifyLessThanOrEqual(exactsolu1(1),result.intervals(1));
            testCase.verifyLessThanOrEqual(result.intervals(1),result.intervals(2));
            testCase.verifyLessThanOrEqual(result.intervals(2),exactsolu1(2));
        end
        
    end
end

