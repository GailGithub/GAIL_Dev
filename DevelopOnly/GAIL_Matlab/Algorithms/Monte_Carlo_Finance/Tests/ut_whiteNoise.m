%% UT_WHITENOISE
% fast unit tests for the whiteNoise class
classdef ut_whiteNoise < matlab.unittest.TestCase

  methods(Test)
             
    function wnNoInput(testCase)
       wn = whiteNoise;
       paths = genPaths(wn,1e6);
       testCase.verifyClass(wn,'whiteNoise');
       testCase.verifyLessThanOrEqual(abs(var(paths)-1/12),0.001);
    end
       
    function wnInputType(testCase)
       param.inputType = 'x';
       wn = whiteNoise(param);
       paths = genPaths(wn,rand(1e6,wn.timeDim.nSteps));
       testCase.verifyClass(wn,'whiteNoise');
       testCase.verifyEqual(wn.inputType,'x');
       testCase.verifyLessThanOrEqual(abs(var(paths)-1/12),0.001);
   end
    
   function wnTimeVector(testCase)
      param.timeDim.timeVector = 1:5;
      wn = whiteNoise(param);
      paths = genPaths(wn,1e6);
      testCase.verifyClass(wn,'whiteNoise');
      testCase.verifyEqual(wn.timeDim.timeVector,1:5);
      testCase.verifyLessThanOrEqual(abs(var(paths)-1/12),0.001);
   end
    
   function wnDim(testCase)
      param.timeDim.dim=2;
      wn = whiteNoise(param);
      testCase.verifyClass(wn,'whiteNoise');
      testCase.verifyEqual(wn.timeDim.dim,2);
   end
    
   function wnSampleKind(testCase)
      param.wnParam.sampleKind='Sobol';
      wn = whiteNoise(param);
      paths = genPaths(wn,1e6);
      testCase.verifyClass(wn,'whiteNoise');
      testCase.verifyEqual(wn.wnParam.sampleKind,'Sobol');
      testCase.verifyLessThanOrEqual(abs(var(paths)-1/12),0.001);
   end
    
   
  end
end
