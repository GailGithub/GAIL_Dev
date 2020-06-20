%UT_FUNMIN_G unit test for funmin_g
classdef ut_funmin_g < matlab.unittest.TestCase
    
    methods(Test)
           
        function funmin_gEXM1(testCase)
            in_param.abstol = 10^(-8); 
            in_param.nmax = 10^5; 
            a = 0.027;
            z = 0.86;
            f = @(x) 0.5/a^2.*(-4*a.^2-(x-z).^2-(x-z-a).*abs(x-z-a)...
                +(x-z+a).*abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);
            [fmin, result] = funmin_g(f,in_param); 
            actualmin = -1;
            actualerr = abs(actualmin-fmin);
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
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        

        function funmin_gEXM2(testCase)
            in_param.abstol = 10^(-8); 
            in_param.nmax = 10^5;  
            a = 0.019;
            z = 0.77;
            f = @(x) 0.5/a^2.*(-4*a.^2-(x-z).^2-(x-z-a).*abs(x-z-a)...
                +(x-z+a).*abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);
            [fmin, result] = funmin_g(f,in_param); 
            actualmin = -1;
            actualerr = abs(actualmin-fmin);
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
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
  
        function funmin_gEXM3(testCase)
            in_param.abstol = 10^(-3); 
            in_param.nmax = 10^5;  
            z = 0.5;
            f = @(x) exp(-1./((x-z).^2)); % flat bottom
            [fmin,result] = funmin_g(f,in_param);
            actualmin = 0;
            actualerr = abs(actualmin-fmin);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.errest,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_gEXM4(testCase)
            z = 1;
            f = @(x) 1./((x-z).^2); % discontinuous function
            testCase.verifyError(@()funmin_g(f),'GAIL:funmin_g:yInf');      
        end
        
        function funmin_gEXM5(testCase)
            in_param.abstol = 1e-6;
            in_param.a = -1; 
            in_param.b = 1;
            in_param.nmax = 10^5; 
            f = @(x) x.^4 .* sin(1./((x==0)+x));
            [fmin, result] = funmin_g(f,in_param); 
            actualmin = sin(-1);
            actualerr = abs(actualmin-fmin);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);   
        end
        
        function funmin_gEXM6(testCase)
            in_param.abstol = 1e-6;
            in_param.a = -1; 
            in_param.b = 1;
            in_param.nmax = 10^5; 
            f = @(x) x.^4 .* sin(1./((x==0)+x)) + 10.*x.^2;
            [fmin, result] = funmin_g(f,in_param); 
            actualmin = 0;
            actualerr = abs(actualmin-fmin);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_gEXM7(testCase)
            in_param.abstol = 1e-6;
            in_param.a = -1; 
            in_param.b = 1;
            in_param.nmax = 10^5; 
            delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
            f = @(x) -B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
                - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta);
            [fmin, result] = funmin_g(f,in_param); 
            actualmin = -1;
            actualerr = abs(actualmin-fmin);
            intnum = size(result.intervals,2);
            exactsolu = c;
            testflag = 1;
            for k=1:intnum
                if exactsolu <= result.intervals(2,k) && exactsolu >= ...
                        result.intervals(1,k)
                 testflag = 0;
                end
            end
            testCase.verifyLessThan(testflag,1);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end
        
        function funmin_gEXM8(testCase)
            in_param.abstol = 1e-6;
            in_param.a = -1; 
            in_param.b = 1;
            in_param.nmax = 10^5; 
            f = @(x) sign(x); 
            warning('off','GAIL:funmin_g:exceediter');
            [fmin, result] = funmin_g(f,in_param); 
            warning('on','GAIL:funmin_g:exceediter');
            actualmin = -1;
            actualerr = abs(actualmin-fmin);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
            testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
        end        
        
    end
end


