%% SpreadOptionTest
% fast unit tests for generating payoffs for spread options
% Author: Tianpei Qian
classdef SpreadOptionTest < matlab.unittest.TestCase

  methods(Test)
             
    function testIIDSolution(testCase)
       
       inp.payoffParam.optType = {'spread'};
       inp.payoffParam.putCallType = {'put'};
       inp.payoffParam.strike = 10;
       inp.assetParam.nAsset = 2;
       inp.assetParam.initPrice = [15 15];
       inp.assetParam.corrMat = [1 1;1 1];
       inp.assetParam.volatility = [0.5 0.5];
       inp.priceParam.absTol = 0.01;
       inp.priceParam.relTol = 0;
       inp.timeDim.timeVector = 1;      
       basket = optPrice(inp);
       price = genOptPrice(basket);
       testCase.verifyLessThanOrEqual(abs(price-10*exp(-0.01)),0.01);
    end
    
        function testSobolSolution(testCase)
       
       inp.payoffParam.optType = {'spread'};
       inp.payoffParam.putCallType = {'put'};
       inp.payoffParam.strike = 10;
       inp.assetParam.nAsset = 2;
       inp.assetParam.initPrice = [15 15];
       inp.assetParam.corrMat = [1 1;1 1];
       inp.assetParam.volatility = [0.5 0.5];
       inp.priceParam.cubMethod='Sobol';
       inp.priceParam.absTol = 0.01;
       inp.priceParam.relTol = 0;
       inp.timeDim.timeVector = 1;      
       basket = optPrice(inp);
       price = genOptPrice(basket);
       testCase.verifyLessThanOrEqual(abs(price-10*exp(-0.01)),0.01);
    end
  end
end
