%% ut_assetPath
% fast unit tests for the assetPath class
classdef ut_assetPath < matlab.unittest.TestCase

  methods(Test)
             
    function assetNoInput(testCase)
       asset = assetPath;
       paths = genPaths(asset,1e6);
       testCase.verifyClass(asset,'assetPath');
       testCase.verifyLessThanOrEqual( ...
          abs(mean(paths,1) .* exp(-asset.assetParam.interest.*asset.timeDim.timeVector) ...
          ./ asset.assetParam.initPrice - 1),0.01);
    end
       
    function assetInputType(testCase)
       param.inputType = 'x';
       asset = assetPath(param);
       paths = genPaths(asset,rand(1e6,asset.timeDim.nSteps));
       testCase.verifyClass(asset,'assetPath');
       testCase.verifyEqual(asset.inputType,'x');
       testCase.verifyLessThanOrEqual( ...
          abs(mean(paths,1) .* exp(-asset.assetParam.interest.*asset.timeDim.timeVector) ...
          ./ asset.assetParam.initPrice - 1),0.01);
    end
    
    function assetTimeVector(testCase)
       param.timeDim.timeVector = 1:5;
       asset = assetPath(param);
       paths = genPaths(asset,1e6);
       testCase.verifyClass(asset,'assetPath');
       testCase.verifyEqual(asset.timeDim.timeVector,1:5);
       testCase.verifyLessThanOrEqual( ...
          abs(mean(paths,1) .* exp(-asset.assetParam.interest.*asset.timeDim.timeVector) ...
          ./ asset.assetParam.initPrice - 1),0.01);
    end
    
%     function assetDim(testCase)
%        param.timeDim.dim=2;
%        asset = assetPath(param);
%        testCase.verifyClass(asset,'assetPath');
%        testCase.verifyEqual(asset.timeDim.dim,2);
%     end
    
    function assetSampleKind(testCase)
       param.wnParam.sampleKind='Sobol';
       asset = assetPath(param);
       paths = genPaths(asset,1e6);
       testCase.verifyClass(asset,'assetPath');
       testCase.verifyEqual(asset.wnParam.sampleKind,'Sobol');
       testCase.verifyLessThanOrEqual( ...
          abs(mean(paths,1) .* exp(-asset.assetParam.interest.*asset.timeDim.timeVector) ...
          ./ asset.assetParam.initPrice - 1),0.01);
    end
    
   
  end
end
