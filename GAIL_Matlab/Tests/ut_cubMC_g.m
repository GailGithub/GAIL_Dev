%UT_CUBMC_G  unit test for cubMC_g
classdef ut_cubMC_g < matlab.unittest.TestCase
    
    methods(Test)
        
        function cubMC_gOfwarning(testCase)
            testCase.verifyWarning(@()cubMC_g,'GAIL:cubMC_g:fnotgiven');
        end
        function cubMC_gOferror10(testCase)
            testCase.verifyError(@()cubMC_g(@(x)x.^2,nan),...
                'GAIL:cubMC_g:hyperboxnotnum');
        end      
        function cubMC_gOferror11(testCase)
            testCase.verifyError(@()cubMC_g(@(x)x.^2,1),...
                'GAIL:cubMC_g:hyperboxnot2d');
        end
        function cubMC_gOferror12(testCase)
            testCase.verifyError(@()cubMC_g(@(x)x.^2,[1 1]),...
                'GAIL:cubMC_g:hyperboxnotlessthan2');
        end
        function cubMC_gOferror13(testCase)
            testCase.verifyError(@()cubMC_g(@(x)x.^2,[-inf,1],...
                'measure','uniform'),'GAIL:cubMC_g:hyperboxnotfiniteforuniform');
        end
        function cubMC_gOferror14(testCase)
            testCase.verifyError(@()cubMC_g(@(x)x.^2,[0,1],...
                'measure','normal'),'GAIL:cubMC_g:hyperboxnotinffornormal');
        end
        function cubMC_gOfxsquare(testCase)
            f = @(x) x.^2;
            in_param.abstol = 1e-3;
            in_param.reltol = 0;
            interval=[0;1];
            meanf = cubMC_g(f,interval,in_param);
            exactf = 1/3;
            actualerr = abs(meanf-exactf);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end       
        function cubMC_gOfexp(testCase)
            f = @(x) exp(x);
            in_param.abstol = 0;
            in_param.reltol = 1e-2;
            interval=[0;1];
            meanf = cubMC_g(f,interval,in_param);
            exactf = exp(1)-1;
            actualerr = abs(meanf-exactf)/exactf;
            testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);
        end        
        function cubMC_gOfsin(testCase)
            f = @(x) sin(x);
            in_param.abstol = 1e-3;
            in_param.reltol = 1e-13;
            interval=[0;1];
            meanf = cubMC_g(f,interval,in_param);
            exactf = 1-cos(1);
            actualerr = abs(meanf-exactf);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end        
        function cubMC_gOfmultierrfun(testCase)
            f = @(x) exp(-x(:,1).^2-x(:,2).^2);
            in_param.abstol = 1e-3;
            in_param.reltol = 1e-13;
            interval=[0 0;1 1];
            meanf = cubMC_g(f,interval,in_param);
            exactf = pi/4*erf(1)^2;
            actualerr = abs(meanf-exactf);
            testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
        end
        
        function cubMC_gOfWorkouts(testCase)
            [ut_abserr,ut_relerr,abstol,reltol] = Test_cubMC_g;
            verifyabserr = ut_abserr<=abstol;
            verifyrelerr = ut_relerr<=reltol;
            testCase.verifyTrue(min(min(verifyabserr + verifyrelerr))>0);
        end
        
        function cubMC_gNormal(testCase)
            format compact
            warning('off','GAIL:meanMC_g:maxreached')
            f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
            count = 0;
            for i=1:20
                [q,out_param] = cubMC_g(f,hyperbox,'normal',1e-3,1e-3);
                exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
                if check==0 || isfinite(q) ==0
                    i, exactsol, q, exitflag = out_param.exitflag,
                    abserr = abs(exactsol-q), tol = gail.tolfun(1e-3,1e-3,1,exactsol,'max')
                    disp('-----');
                    count = count + 1;
                    %keyboard
                else
                    gail.print_iterations(i,'i',true);
                end
            end
            warning('on','GAIL:meanMC_g:maxreached')
            testCase.verifyTrue(count==0);
        end
    end
end

