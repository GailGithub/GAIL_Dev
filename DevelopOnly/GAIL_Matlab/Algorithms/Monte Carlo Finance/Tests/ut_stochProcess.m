%% ut_stochProcess
% fast unit tests for the stochProcess class
classdef ut_stochProcess < matlab.unittest.TestCase

  methods(Test)
             
    function spNoInput(testCase)
       sp=stochProcess;
       testCase.verifyClass(sp,'stochProcess');
    end
       
    function spInputType(testCase)
       param.inputType='x';
       sp=stochProcess(param);
       testCase.verifyClass(sp,'stochProcess');
       testCase.verifyEqual(sp.inputType,'x');
    end
    
    function spTimeVector(testCase)
       param.timeDim.timeVector=1:5;
       sp=stochProcess(param);
       testCase.verifyClass(sp,'stochProcess');
       testCase.verifyEqual(sp.timeDim.timeVector,1:5);
    end
    
   function spDim(testCase)
      param.timeDim.dim=2;
      sp=stochProcess(param);
      testCase.verifyClass(sp,'stochProcess');
      testCase.verifyEqual(sp.timeDim.dim,2);
   end
    
   
  end
end
