%UT_MEANMC_G  unit tests for meanMC_g and Test_meanMC_g (workouts)

classdef ut_meanMC_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function meanMC_gOfexp(testCase)
      in_param.abstol = 1e-2;
      in_param.reltol = 0;
      meanY = meanMC_g(@(n) exp(rand(n,1)),in_param);
      exactY = exp(1)-1;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
    function meanMC_gOfxsquare(testCase)
      in_param.reltol = 1e-1;
      in_param.abstol = 0;
      meanY = meanMC_g(@(n) rand(n,1).^2,in_param);
      exactY = 1/3;
      actualerr = abs(meanY-exactY)/exactY;
      testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);
    end
    
    function meanMC_gOfsin(testCase)
      in_param.reltol = 1e-2;
      in_param.abstol = 0;
      meanY = meanMC_g(@(n) sin(rand(n,1)),in_param);
      exactY = 1-cos(1);
      actualerr = abs(meanY-exactY)/exactY;
      testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);
    end
    
    function meanMC_gOfparsing(testCase)
      in_param.abstol = -1e-2;  
      in_param.reltol = 0;
      meanY = testCase.verifyWarning(@()meanMC_g(@(n) rand(n,1).^2,...
        in_param),'GAIL:meanMC_g:abstolneg');
      exactY = 1/3;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,abs(in_param.abstol));
    end
    
    function meanMC_gOfnonRandomInput(testCase)
        in_param.abstol = 1e-2;
        testCase.verifyWarning(@()meanMC_g(@(x) x.^2,...
            in_param),'GAIL:meanMC_g:yrandnotlengthN');  
    end
    function meanMC_gOfWorkouts(testCase)
         mu = Test_meanMC_g;
        testCase.verifyTrue(mu>4.4);
        testCase.verifyTrue(mu<4.5);  
    end
  end
end

