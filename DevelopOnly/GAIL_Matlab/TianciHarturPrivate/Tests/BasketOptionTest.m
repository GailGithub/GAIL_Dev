  classdef BasketOptionTest < matlab.unittest.TestCase
      % BasketOptionTest tests the basket option price generated using
      % Monte Carlo Method and Quasi Monte Carlo Method.
    
      methods (Test)
          function testIIDSolution(testCase)
              inp.payoffParam.optType = {'basket'};
              %inp.payoffParam.basketWeight = [0.5 0.5];
              %inp.assetParam.initPrice = [11 11]; 
              %inp.assetParam.volatility = [0.5 0.5];
              %inp.assetParam.nAsset = 2;
              %inp.assetParam.sqCorr = [1 1;1 1];
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              BasketOption = optPrice(inp);
              [expSolutionCall,expSolutionPut] = vanillaPrice(BasketOption.assetParam.initPrice,...
                  BasketOption.payoffParam.strike, BasketOption.assetParam.interest,...
                  BasketOption.assetParam.volatility, BasketOption.timeDim.endTime,...
                  BasketOption.payoffParam.basketWeight);
              actSolutionCall = genOptPrice(BasketOption);
              BasketOption.payoffParam.putCallType = {'put'};
              actSolutionPut = genOptPrice(BasketOption);
              testCase.verifyLessThan(abs(actSolutionCall-expSolutionCall),0.01);
              testCase.verifyLessThan(abs(actSolutionPut-expSolutionPut),0.01);
          end
          function testSobolSolution(testCase)
              inp.payoffParam.optType = {'basket'};
              %inp.payoffParam.basketWeight = [0.5 0.5];
              %inp.assetParam.initPrice = [11 11]; 
              %inp.assetParam.volatility = [0.5 0.5];
              %inp.assetParam.nAsset = 2;
              %inp.assetParam.sqCorr = [1 1;1 1];
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              inp.priceParam.cubMethod = 'Sobol';
              BasketOption = optPrice(inp);
              [expSolutionCall,expSolutionPut] = vanillaPrice(BasketOption.assetParam.initPrice,...
                  BasketOption.payoffParam.strike, BasketOption.assetParam.interest,...
                  BasketOption.assetParam.volatility, BasketOption.timeDim.endTime,...
                  BasketOption.payoffParam.basketWeight);
              actSolutionCall = genOptPrice(BasketOption);
              BasketOption.payoffParam.putCallType = {'put'};
              actSolutionPut = genOptPrice(BasketOption);
              testCase.verifyLessThan(abs(actSolutionCall-expSolutionCall),0.01);
              testCase.verifyLessThan(abs(actSolutionPut-expSolutionPut),0.01);
          end
      end
      
  end 