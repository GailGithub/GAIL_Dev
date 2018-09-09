%UT_CUBBAYESLATTICE_G  unit test for cubBayesLattice_g
classdef ut_cubBayesLattice_g < matlab.unittest.TestCase
    
    methods(Test)
        
        function cubBayesLattice_gOfwarning(testCase)
            testCase.verifyWarning(@()cubBayesLattice_g(),...
              'GAIL:cubBayesLattice_g:fdnotgiven');
        end
        function cubBayesLattice_gOferror10(testCase)
            testCase.verifyWarning(@()cubBayesLattice_g('f',@(x)x.^2),...
                'GAIL:cubBayesLattice_g:fdnotgiven');
        end
        function cubBayesLattice_gOferror11(testCase)
            testCase.verifyWarning(@()cubBayesLattice_g('f',@(x)x.^2,...
              'dim',2,'order',4),...
                'GAIL:cubBayesLattice_g:r_invalid');
        end
        function cubBayesLattice_gOferror12(testCase)
            testCase.verifyWarning(@()cubBayesLattice_g('f',@(x)x.^2,'dim',2,...
              'stopCriterion','XFZ'),...
                'GAIL:cubBayesLattice_g:stop_crit_invalid');
        end
        function cubBayesLattice_gOferror13(testCase)
            testCase.verifyWarning(@()cubBayesLattice_g('f',@(x)x.^2,'dim',2,...
                'ptransform','uniform'),...
                'GAIL:cubBayesLattice_g:var_transform_invalid');
        end

        function cubBayesLattice_gOfxsquare(testCase)
            f = @(x) x.^2;
            abstol = 1e-3;
            reltol = 0;
            dim=1;
            obj = cubBayesLattice_g('f',f,'dim',dim,'absTol',abstol,'reltol',reltol);
            meanf = compInteg(obj);
            exactf = 1/3;
            actualerr = abs(meanf-exactf);
            testCase.verifyLessThanOrEqual(actualerr,abstol);
        end       
        function cubBayesLattice_gOfexp(testCase)
            f = @(x) exp(x);
            abstol = 0;
            reltol = 1e-2;
            dim=1;
            obj = cubBayesLattice_g('f',f,'dim',dim,'absTol',abstol,...
              'reltol',reltol); %, 'order',2, 'ptransform','C1sin'
            meanf = compInteg(obj);
            exactf = exp(1)-1;
            actualerr = abs(meanf-exactf)/exactf;
            testCase.verifyLessThanOrEqual(actualerr,reltol);
        end        
        function cubBayesLattice_gOfsin(testCase)
            f = @(x) sin(x);
            abstol = 1e-3;
            reltol = 1e-13;
            dim=1;
            obj = cubBayesLattice_g('f',f,'dim',dim,'absTol',abstol,...
              'reltol',reltol);
            meanf = compInteg(obj);
            exactf = 1-cos(1);
            actualerr = abs(meanf-exactf);
            testCase.verifyLessThanOrEqual(actualerr,abstol);
        end        
        function cubBayesLattice_gOfmultierrfun(testCase)
            f = @(x) exp(-x(:,1).^2-x(:,2).^2);
            abstol = 1e-3;
            reltol = 1e-13;
            dim=2;
            obj = cubBayesLattice_g('f',f,'dim',dim,'absTol',abstol,...
              'reltol',reltol);
            meanf = compInteg(obj);
            exactf = pi/4*erf(1)^2;
            actualerr = abs(meanf-exactf);
            testCase.verifyLessThanOrEqual(actualerr,abstol);
        end
        
%         function cubBayesLattice_gOfWorkouts(testCase)
%             [ut_abserr,ut_relerr,abstol,reltol] = Test_cubBayesLattice_g;
%             verifyabserr = ut_abserr<=abstol;
%             verifyrelerr = ut_relerr<=reltol;
%             testCase.verifyTrue(min(min(verifyabserr + verifyrelerr))>0);
%         end
%         
%         function cubBayesLattice_gNormal(testCase)
%             format compact
%             warning('off','GAIL:meanMC_g:maxreached')
%             f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];dim=3;
%             count = 0;
%             for i=1:20
%                 [q,out_param] = cubBayesLattice_g(f,hyperbox,'normal',1e-3,1e-3);
%                 exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
%                 if check==0 || isfinite(q) ==0,
%                     i, exactsol, q, exitflag = out_param.exitflag,
%                     abserr = abs(exactsol-q), tol = gail.tolfun(1e-3,1e-3,1,exactsol,'max')
%                     disp('-----');
%                     count = count + 1;
%                     %keyboard
%                 else
%                     i
%                 end;
%             end;
%             warning('on','GAIL:meanMC_g:maxreached')
%             testCase.verifyTrue(count==0);
%         end
    end
end

