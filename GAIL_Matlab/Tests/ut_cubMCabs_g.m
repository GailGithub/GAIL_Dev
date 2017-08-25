%UT_CUBMCabs_G  unit test for cubMCabs_g
classdef ut_cubMCabs_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubMCabs_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 1e-2;
      interval=[0;1];
      meanf = cubMCabs_g(f,interval,in_param);
      exactf = 1/3;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubMCabs_gOfexp(testCase)
      f = @(x) exp(x);
      in_param.abstol = 1e-3;
      interval=[0;1];
      meanf = cubMCabs_g(f,interval,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubMCabs_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 1e-3;
      interval=[0;1];
      meanf = cubMCabs_g(f,interval,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubMCabs_gOfmultierrfun(testCase)
      f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      in_param.abstol = 1e-3;
      interval=[0 0;1 1];
      meanf = cubMCabs_g(f,interval,in_param);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function cubMCabs_gOfwarning(testCase)
        testCase.verifyWarning(@()cubMCabs_g,'GAIL:cubMCabs_g:fnotgiven');
    end
    function cubMCabs_gOferror10(testCase)
        testCase.verifyError(@()cubMCabs_g(@(x)x.^2,nan),...
            'GAIL:cubMCabs_g:hyperboxnotnum');
    end
    
    function cubMCabs_gOferror11(testCase)
        testCase.verifyError(@()cubMCabs_g(@(x)x.^2,1),...
            'GAIL:cubMCabs_g:hyperboxnot2d');
    end
    function cubMCabs_gOferror12(testCase)
        testCase.verifyError(@()cubMCabs_g(@(x)x.^2,[1 1]),...
            'GAIL:cubMCabs_g:hyperboxnotlessthan2');
    end
    function cubMCabs_gOferror13(testCase)
        testCase.verifyError(@()cubMCabs_g(@(x)x.^2,[-inf,1],...
            'measure','uniform'),'GAIL:cubMCabs_g:hyperboxnotfiniteforuniform');
    end
    function cubMCabs_gOferror14(testCase)
        testCase.verifyError(@()cubMCabs_g(@(x)x.^2,[0,1],...
            'measure','normal'),'GAIL:cubMCabs_g:hyperboxnotinffornormal');
    end
  end
end

