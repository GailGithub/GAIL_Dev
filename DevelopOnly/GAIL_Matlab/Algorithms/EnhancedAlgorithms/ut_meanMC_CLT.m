%UT_meanMC_CLT  unit tests for meanMC_CLT 

classdef ut_meanMC_CLT < matlab.unittest.TestCase
  
  methods(Test)
    
    function meanMC_CLTOfexp(testCase)
      in_param.abstol = 1e-2;
      in_param.reltol = 0;
      meanY = meanMC_CLT(@(n) exp(rand(n,1)),in_param.abstol,in_param.reltol);
      exactY = exp(1)-1;
      actualerr = abs(meanY-exactY);
      testCase.verifyLessThanOrEqual(actualerr,in_param.abstol);
    end
    
     function meanMC_CLTOfxsquare(testCase)

      in_param.reltol = 1e-1;

      in_param.abstol = 0;

      meanY = meanMC_CLT(@(n) rand(n,1).^2,in_param.abstol,in_param.reltol);

      exactY = 1/3;

      actualerr = abs(meanY-exactY)/exactY;

      testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);

     end
    
     function meanMC_CLTOfsin(testCase)

      in_param.reltol = 1e-2;

      in_param.abstol = 0;

      meanY = meanMC_CLT(@(n) sin(rand(n,1)),in_param.abstol,in_param.reltol);

      exactY = 1-cos(1);

      actualerr = abs(meanY-exactY)/exactY;

      testCase.verifyLessThanOrEqual(actualerr,in_param.reltol);

    end

    function meanMC_CLTOfparsing(testCase)
      in_param.abstol = -1e-2;  

      in_param.reltol = 0;
      errorOccurred = false;
      
      try
         meanMC_CLT(@(n) rand(n,1).^2,in_param.abstol,in_param.reltol)
      catch ME
         testCase.verifyEqual(ME.identifier,'MATLAB:expectedNonnegative')
         errorOccurred = true;
      end
      testCase.verifyEqual(errorOccurred,true);


    end



  function meanMC_CLTOfnonRandomInput(testCase)

        in_param.abstol = 1e-2;
        errorOccurred = false;
        try
           a=meanMC_CLT(@(x) x.^2,in_param.abstol);
        catch ME
           display('a')
           testCase.verifyEqual(ME.identifier,'MATLAB:expectedNonnegative')
           errorOccurred = true;
        end
        testCase.verifyEqual(errorOccurred,true);
        
   end

    %function meanMC_CLTOfWorkouts(testCase)

   %      mu = Test_meanMC_CLT;

    %   testCase.verifyTrue(mu>4.4);

     %   testCase.verifyTrue(mu<4.5);  

  %  end



   
  end
end
