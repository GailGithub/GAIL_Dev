%ut_cubLattice_g  unit test for cubLattice_g
classdef ut_cubLattice_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubLattice_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 1e-2;
      d=1;
      meanf = cubLattice_g(f,d,in_param);
      exactf = 1+0.33;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubLattice_gOfexp(testCase)
      f = @(x) exp(x);
      in_param.abstol = 1e-3;
      d=1;
      meanf = cubLattice_g(f,d,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubLattice_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 1e-3;
      d=1;
      meanf = cubLattice_g(f,d,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubLattice_gOfmultierrfun(testCase)
      f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      in_param.abstol = 1e-3;
      d=2;
      meanf = cubLattice_g(f,d,in_param);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubLattice_gOfwarning(testCase)
        testCase.verifyWarning(@()cubLattice_g,'MATLAB:cubLattice_g:fdnotgiven');
    end
    function cubLattice_gOdwarning(testCase)
        testCase.verifyWarning(@()cubLattice_g(@(x)x.^2,1.5),'MATLAB:cubLattice_g:dnotposint');
    end
  end
end
