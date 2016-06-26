%% OutperformanceOptionTest
% fast unit tests for generating payoffs for outperformance options
% Author: Tianpei Qian
classdef OutperformanceOptionTest < matlab.unittest.TestCase

  methods(Test)
             
    function testIIDSolution(testCase)
       addpath('../OptionPricing') % folder of optPayoff
       
       % parameter for the euro call option
       euroinp.payoffParam.strike = 10;
       euroinp.assetParam.initPrice = 15;
       euroinp.assetParam.volatility = 0.5;
       euroinp.timeDim.timeVector = 1;      
       euro = optPayoff(euroinp);
       exactPrice = euro.exactPrice;
       
       
       inp.payoffParam.optType = {'outperformance'};
       inp.payoffParam.strike = 10;
       inp.payoffParam.outperformanceWeight = [1 1];
       inp.assetParam.nAsset = 2;
       inp.assetParam.initPrice = [15 15];
       inp.assetParam.corrMat = [1 1;1 1];
       inp.assetParam.volatility = [0.5 0.5];
       inp.priceParam.absTol = 0.01;
       inp.priceParam.relTol = 0;
       inp.timeDim.timeVector = 1;      
       outperform = optPrice(inp);
       price = genOptPrice(outperform);
       testCase.verifyLessThanOrEqual(abs(price-exactPrice),0.01);
    end
    
    function testSobolSolution(testCase)
       addpath('../OptionPricing') % folder of optPayoff
       
       % parameter for the euro call option
       euroinp.payoffParam.strike = 10;
       euroinp.assetParam.initPrice = 15;
       euroinp.assetParam.volatility = 0.5;
       euroinp.timeDim.timeVector = 1;      
       euro = optPayoff(euroinp);
       exactPrice = euro.exactPrice;
       
       
       inp.payoffParam.optType = {'outperformance'};
       inp.payoffParam.strike = 10;
       inp.payoffParam.outperformanceWeight = [1 1];
       inp.assetParam.nAsset = 2;
       inp.assetParam.initPrice = [15 15];
       inp.assetParam.corrMat = [1 1;1 1];
       inp.assetParam.volatility = [0.5 0.5];
       inp.priceParam.cubMethod='Sobol';
       inp.priceParam.absTol = 0.01;
       inp.priceParam.relTol = 0;
       inp.timeDim.timeVector = 1;      
       outperform = optPrice(inp);
       price = genOptPrice(outperform);
       testCase.verifyLessThanOrEqual(abs(price-exactPrice),0.01);
    end
    
  end
end
