%ut_cubLattice_g  unit test for cubLattice_g
classdef ut_cubMC_CLT < matlab.unittest.TestCase
   
   methods(Test)
      
      function cubLattice_gOfxsquare(testCase)
         w.f= @(x) x.^2;
         w.absTol=1e-2;
         w.domain = [0;1];
         [meanf, out_param]=cubMC_CLT(w);
         exactf = 0.33;
         actualerr = abs(meanf-exactf);
         tolerance = max(out_param.err.absTol,out_param.err.relTol*abs(exactf));
         testCase.verifyLessThanOrEqual(actualerr,tolerance);
      end
      
      function cubLattice_gOfexp(testCase)
         w.f = @(x) exp(x);
         w.absTol = 1e-3;
         w.domain = [0;1];
         [meanf, out_param]=cubMC_CLT(w);
         exactf = exp(1)-1;
         actualerr = abs(meanf-exactf);
         tolerance = max(out_param.err.absTol,out_param.err.relTol*abs(exactf));
         testCase.verifyLessThanOrEqual(actualerr,tolerance);
         testCase.verifyTrue(out_param.fun.d==1);
      end
      
      function cubLattice_gOfsin(testCase)
         w.f = @(x) sin(x);
         w.absTol = 1e-3;
         w.domain = [0;1];
         [meanf, out_param]=cubMC_CLT(w);
         exactf = 1-cos(1);
         actualerr = abs(meanf-exactf);
         tolerance = max(out_param.err.absTol,out_param.err.relTol*abs(exactf));
         testCase.verifyLessThanOrEqual(actualerr,tolerance);
         testCase.verifyTrue(out_param.fun.d==1);
      end
      
      function cubLattice_gOfmultierrfun(testCase)
         w.f = @(x) exp(-x(:,1).^2-x(:,2).^2);
         w.absTol = 1e-3;
         w.domain = [0 0;1 1];
         [meanf, out_param]=cubMC_CLT(w);
         exactf = pi/4*erf(1)^2;
         actualerr = abs(meanf-exactf);
         tolerance = max(out_param.err.absTol,out_param.err.relTol*abs(exactf));
         testCase.verifyLessThanOrEqual(actualerr,tolerance);
         testCase.verifyTrue(out_param.fun.d==2);
      end
      
      function cubLattice_gOfwarning(testCase)
         testCase.verifyWarning(@()cubLattice_g,'GAIL:cubLattice_g:fdnotgiven');
      end
      
      function cubLattice_gOdwarning(testCase)
         testCase.verifyWarning(@()cubLattice_g(@(x)x.^2,1.5),'GAIL:cubLattice_g:hyperbox_error1');
      end
            
      function cubLattice_Workouts(testCase)
         [ut_abserr,ut_relerr,abstol,reltol] = Test_cubLattice_g;
         verifyabserr = ut_abserr<=abstol;
         verifyrelerr = ut_relerr<=reltol;
         testCase.verifyTrue(min(min(verifyabserr + verifyrelerr))>0);
      end
   end
end
