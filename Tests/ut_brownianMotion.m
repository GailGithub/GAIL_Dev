%% ut_brownianMotion
% fast unit tests for the brownianMotion class
classdef ut_brownianMotion < matlab.unittest.TestCase

  methods(Test)
             
    function bmNoInput(testCase)
       bm = brownianMotion;
       paths = genPaths(bm,1e6);
       testCase.verifyClass(bm,'brownianMotion');
       testCase.verifyLessThanOrEqual( ...
          abs(var(paths)-bm.timeDim.timeVector)./bm.timeDim.timeVector,0.01);
    end
       
    function bmInputType(testCase)
       param.inputType = 'x';
       bm = brownianMotion(param);
       paths = genPaths(bm,rand(1e6,bm.timeDim.nSteps));
       testCase.verifyClass(bm,'brownianMotion');
       testCase.verifyEqual(bm.inputType,'x');
       testCase.verifyLessThanOrEqual( ...
          abs(var(paths)-bm.timeDim.timeVector)./bm.timeDim.timeVector,0.01);
   end
    
   function bmTimeVector(testCase)
      param.timeDim.timeVector = 1:5;
      bm = brownianMotion(param);
      paths = genPaths(bm,1e6);
      testCase.verifyClass(bm,'brownianMotion');
      testCase.verifyEqual(bm.timeDim.timeVector,1:5);
      testCase.verifyLessThanOrEqual( ...
         abs(var(paths)-bm.timeDim.timeVector)./bm.timeDim.timeVector,0.01);
   end
    
   function bmDim(testCase)
      param.bmParam.whBM = [1 2];
      bm = brownianMotion(param);
      testCase.verifyClass(bm,'brownianMotion');
      testCase.verifyEqual(bm.timeDim.dim,2);
   end
    
   function bmSampleKind(testCase)
      param.wnParam.sampleKind='Sobol';
      bm = brownianMotion(param);
      paths = genPaths(bm,1e6);
      testCase.verifyClass(bm,'brownianMotion');
      testCase.verifyEqual(bm.wnParam.sampleKind,'Sobol');
      testCase.verifyLessThanOrEqual( ...
         abs(var(paths)-bm.timeDim.timeVector),0.01);
   end
    
   
  end
end


