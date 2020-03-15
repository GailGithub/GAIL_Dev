%UT_CUBBAYESNET_G  unit tests for cubBayesNet_g
classdef ut_cubBayesNet_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubBayesNet_gOfwarning(testCase)
      testCase.verifyWarning(@()cubBayesNet_g(),...
        'GAIL:cubBayesNet_g:fdnotgiven');
    end
    
    function cubBayesNet_gOferror10(testCase)
      testCase.verifyWarning(@()cubBayesNet_g('f',@(x)x.^2),...
        'GAIL:cubBayesNet_g:fdnotgiven');
    end
    
    function cubBayesNet_gOferror11(testCase)
      testCase.verifyWarning(@()cubBayesNet_g('f',@(x)x.^2,...
        'dim',2,'order',4),...
        'GAIL:cubBayesNet_g:r_invalid');
    end
    
    function cubBayesNet_gOferror12(testCase)
      testCase.verifyWarning(@()cubBayesNet_g('f',@(x)x.^2,'dim',2,...
        'stopCriterion','XFZ'),...
        'GAIL:cubBayesNet_g:stop_crit_invalid');
    end
    
    function cubBayesNet_gOfwarningMaxReached(testCase)
      testCase.verifyWarning( ...
      @()TestAsianArithmeticMeanOptionAutoExample('absTol',1E-7, ...
        'order',1, 'stopAtTol',true, ....
        'stopCriterion','GCV', 'samplingMethod','Net', ...
        'nRepAuto',1, 'log10ErrVec', -7:1:-4 ), ...
        'GAIL:cubBayesNet_g:maxreached');
    end
    
    function cubBayesNet_gOfxsquare(testCase)
      f = @(x) x.^2;
      abstol = 1e-3;
      reltol = 0;
      dim=1;
      obj = cubBayesNet_g('f',f,'dim',dim,'absTol',abstol,...
        'relTol',reltol);
      meanf = compInteg(obj);
      exactf = 1/3;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,abstol);
    end
    
    function cubBayesNet_gOfexp(testCase)
      f = @(x) exp(x);
      abstol = 0;
      reltol = 1e-2;
      dim=1;
      obj = cubBayesNet_g('f',f,'dim',dim,'absTol',abstol,...
        'relTol',reltol);
      meanf = compInteg(obj);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf)/exactf;
      testCase.verifyLessThanOrEqual(actualerr,reltol);
    end
    
    function cubBayesNet_gOfsin(testCase)
      f = @(x) sin(x);
      abstol = 1e-3;
      reltol = 1e-13;
      dim=1;
      obj = cubBayesNet_g('f',f,'dim',dim,'absTol',abstol,...
        'relTol',reltol);
      meanf = compInteg(obj);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,abstol);
    end
    
    function cubBayesNet_gOfmultierrfun(testCase)
      f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      abstol = 1e-3;
      reltol = 1e-13;
      dim=2;
      obj = cubBayesNet_g('f',f,'dim',dim,'absTol',abstol,...
        'relTol',reltol);
      meanf = compInteg(obj);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      testCase.verifyLessThanOrEqual(actualerr,abstol);
    end
    
    function cubBayesNet_gOptPrice(testCase)
      absTolArg = 1E-5;
      inputArgs = {'absTol',absTolArg, 'order',1, ....
        'stopAtTol',true, 'stopCriterion','GCV', ...
        'samplingMethod','Net', 'nRepAuto',1, ...
        'log10ErrVec', -7:1:-4};
      
      format compact
      warning('off','GAIL:cubBayesNet_g:maxreached')
      [muhat,actualerr] = ...
        TestAsianArithmeticMeanOptionAutoExample(inputArgs{:});
      warning('on','GAIL:cubBayesNet_g:maxreached')
      testCase.verifyWarning(@()cubBayesNet_g(),...
        'GAIL:cubBayesNet_g:fdnotgiven');
      testCase.verifyLessThanOrEqual(absTolArg,actualerr);
    end
    
  end
end

