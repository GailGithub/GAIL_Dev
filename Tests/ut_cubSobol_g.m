%ut_cubSobol_g  unit test for cubSobol_g
classdef ut_cubSobol_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubSobol_gOfxsquare(testCase)
      f = @(x) x.^2;
      in_param.abstol = 1e-2;
      hyperbox = [0;1];
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      exactf = 0.33;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubSobol_gOfexp(testCase)
      f = @(x) exp(x);
      in_param.abstol = 1e-3;
      hyperbox = [0;1];
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      exactf = exp(1)-1;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubSobol_gOfsin(testCase)
      f = @(x) sin(x);
      in_param.abstol = 1e-3;
      hyperbox = [0;1];
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      exactf = 1-cos(1);
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==1);
    end
    
    function cubSobol_gOfmultierrfun(testCase)
      f = @(x) exp(-x(:,1).^2-x(:,2).^2);
      in_param.abstol = 1e-3;
      hyperbox = [0 0;1 1];
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==2);
    end

    function cubSobol_gOfpolycv(testCase)
      f.func = @(x) [10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
      f.cv = [1,4/3]; 
      hyperbox= [zeros(1,3);2*ones(1,3)];
      in_param.abstol = 1e-4;
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      exactf = 16/3;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==3);
    end
    
    function cubSobol_gOfmultierfcv(testCase)
      f.func=@(x) [exp(-x(:,1).^2-x(:,2).^2), exp(-x(:,2).^2), x(:,1).^2.*exp(-x(:,2).^2), x(:,1).^4.*exp(-x(:,2).^2)];
      f.cv=[0.5*sqrt(pi)*erf(1), 1/6*sqrt(pi)*erf(1), 0.1*sqrt(pi)*erf(1)]; 
      in_param.abstol = 1e-4;
      hyperbox = [0 0;1 1];
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      exactf = pi/4*erf(1)^2;
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==2);
    end

    function cubSobol_gOfgmeanOptcv(testCase)
      % set up option params
      inp.timeDim.timeVector = 1/4:1/4:4/4;%weekly monitor for a month 
      inp.assetParam.initPrice = 120; %initial stock price
      inp.assetParam.interest = 0.01; %risk-free interest rate
      inp.assetParam.volatility = 0.5; %volatility
      inp.payoffParam.strike = 130; %strike price
      inp.payoffParam.barrier = 135; %strike price
      inp.priceParam.absTol = 1e-3; %absolute tolerance of a nickel
      inp.priceParam.relTol = 0; %zero relative tolerance
      inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
      inp.bmParam.assembleType = 'diff';
      obj = optPrice(inp); %construct an optPrice object
      opt = optPayoff(obj);
      opt.payoffParam = struct( ...
      'optType',{{'gmean','euro'}},...
      'putCallType',{{'call','call'}}); 
      exactOpt = opt.exactPrice;    
      % define the function
      f.func = @(x) genOptPayoffs(opt,x);
      f.cv = exactOpt(2:end); 
      in_param.abstol = inp.priceParam.absTol;
      in_param.reltol = inp.priceParam.relTol;
      hyperbox = [zeros(1,opt.timeDim.nSteps);ones(1,opt.timeDim.nSteps)];
      exactf = exactOpt(1);
      [meanf,out_param] = cubSobol_g(f,hyperbox,in_param);
      actualerr = abs(meanf-exactf);
      tolerance = max(out_param.abstol,out_param.reltol*abs(exactf));
      testCase.verifyLessThanOrEqual(actualerr,tolerance);
      testCase.verifyTrue(out_param.d==4);
    end

    function cubSobol_gOfwarning(testCase)
        testCase.verifyWarning(@()cubSobol_g,'GAIL:cubSobol_g:fdnotgiven');
    end
    
    function cubSobol_gOdwarning(testCase)
        testCase.verifyWarning(@()cubSobol_g(@(x)x.^2,1.5),'GAIL:cubSobol_g:hyperbox_error1');
    end
    
    function cubSobol_Workouts(testCase)
        [ut_abserr,ut_relerr,abstol,reltol] = Test_cubSobol_g;
        verifyabserr = ut_abserr<=abstol;
        verifyrelerr = ut_relerr<=reltol;
        testCase.verifyTrue(min(min(verifyabserr + verifyrelerr))>0);
    end

  end
end

