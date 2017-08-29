%ut_cubLattice_g  unit test for cubLattice_g
classdef ut_cubLattice_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubLattice_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 1e-2;
      hyperbox = [0;1];
      [meanf,out_param] = cubLattice_g(f,hyperbox,in_param);
      exactf = 0.33;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
    end
    
    function cubLattice_gOfexp(testCase)
      f = @(x) exp(x);
      in_param.abstol = 1e-3;
      hyperbox = [0;1];
      [meanf,out_param] = cubLattice_g(f,hyperbox,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubLattice_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 1e-3;
      hyperbox = [0;1];
      [meanf,out_param] = cubLattice_g(f,hyperbox,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubLattice_gOfmultierrfun(testCase)
      f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      in_param.abstol = 1e-3;
      hyperbox = [0 0;1 1];
      [meanf,out_param] = cubLattice_g(f,hyperbox,in_param);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==2);
    end
    
    function cubLattice_gOfwarning(testCase)
        testCase.verifyWarning(@()cubLattice_g,'GAIL:cubLattice_g:fdnotgiven');
    end
    
    function cubLattice_gOdwarning(testCase)
        testCase.verifyWarning(@()cubLattice_g(@(x)x.^2,1.5),'GAIL:cubLattice_g:hyperbox_error1');
    end
    
    function cubLattice_gNormal(testCase)
        format compact
        f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
        count = 0;
        for i=1:100
            [q,out_param] = cubLattice_g(f,hyperbox,'normal',1e-3,1e-3,'transform','C1sin','shift',2^(-25)*ones(1,3));
            exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
            if check==0 || isfinite(q) ==0 %|| out_param.exitflag > 0,
                i, exactsol, q, exitflag = out_param.exitflag,
                abserr = abs(exactsol-q), tol = gail.tolfun(1e-3,1e-3,1,exactsol,'max'), n = out_param.n,
                shift = out_param.shift, lattice = mod(bsxfun(@plus, gail.lattice_gen(1,2^24,3), shift),1);
                max_lattice = max(max(lattice))
                max_C1sin = max(max(lattice-sin(2*pi*lattice)/(2*pi))),
                max_after_normtransform = max(max(gail.stdnorminv(lattice-sin(2*pi*lattice)/(2*pi))))%, min(min(gail.stdnorminv(lattice-sin(2*pi*lattice)/(2*pi))))
                disp('-----');
                count = count + 1;
                %keyboard
                clear lattice
            end;
        end;
        testCase.verifyTrue(count==0);
    end
    
    function cubLattice_Workouts(testCase)
        [ut_abserr,ut_relerr,abstol,reltol] = Test_cubLattice_g;
        verifyabserr = ut_abserr<=abstol;
        verifyrelerr = ut_relerr<=reltol;
        testCase.verifyTrue(min(min(verifyabserr + verifyrelerr))>0);
    end
  end
end

