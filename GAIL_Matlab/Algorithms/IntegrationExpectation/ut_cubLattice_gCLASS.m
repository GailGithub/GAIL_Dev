%ut_cubLattice_g  unit test for cubLattice_g
classdef ut_cubLattice_gCLASS < matlab.unittest.TestCase
  
  methods(Test)
 
    function cubLattice_gCLASSOfxsquare(testCase)
      w.f = @(x) x.^2;
      w.abstol = 1e-2;
      w.domain = [0;1];
      [meanf,out_param] = cubLattice_gCLASS(w);
      exactf = 0.33;
      actualerr = abs(meanf-exactf);
      
      % Not working 
      tolerance = max(out_param.absTol,out_param.relTol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
    end
    
    function cubLattice_gCLASSOfexp(testCase)
      a.f = @(x) exp(x);
      a.abstol = 1e-3;
      a.domain = [0;1];
      [meanf,out_param] = cubLattice_gCLASS(a);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.absTol,out_param.relTol*abs(exactf));
      
      % Not working
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubLattice_gOfsin(testCase)
      b.f = @(x) sin(x);
      b.abstol = 1e-3;
      b.domain = [0;1];
      [meanf,out_param] = cubLattice_gCLASS(b);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.absTol,out_param.relTol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubLattice_gOfmultierrfun(testCase)
      c.f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      c.abstol = 1e-3;
      c.domain = [0 0;1 1];
      [meanf,out_param] = cubLattice_gCLASS(c);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.absTol,out_param.relTol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==2);
    end
    
    %%%%%%%%%%%%%%%%%%
    function cubLattice_gOfwarning(testCase)
       
       d=cubLattice_gCLASS;
        try
            cubLattice_gCLASS(d);
        catch ME
            testCase.verifyEqual(lastwarn, 'No function input, default used');
            errorOccurred=true;
        end 
            
    end
    
    
    
    function cubLattice_gOdwarning(testCase)
        try
            @()cubLattice_gCLASS(@(x)x.^2,1.5);
        catch ME
            display(ME.identifier);
        end
    end
    
    %%%%%%%%%%%%%%%%%%
    
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
