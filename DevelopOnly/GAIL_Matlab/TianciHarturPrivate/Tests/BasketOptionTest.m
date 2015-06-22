  classdef BasketOptionTest < matlab.unittest.TestCase
      % BasketOptionTest tests the basket option price generated using
      % Monte Carlo Method and Quasi Monte Carlo Method.
    
      methods (Test)
          function testIIDSolution(testCase)
              inp.payoffParam.optType = {'basket'};
              inp.payoffParam.basketWeight = [0.2 0.8];
              inp.assetParam.initPrice = [11 15]; 
              inp.assetParam.volatility = [0.5 0.6];
              inp.assetParam.nAsset = 2;
              inp.assetParam.sqCorr = [1 1; 1 1];
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              BasketOptionCall = optPrice(inp);
              [expSolutionCall,expSolutionPut] = vanillaPrice(BasketOptionCall.assetParam.initPrice,...
                  BasketOptionCall.payoffParam.strike, BasketOptionCall.assetParam.interest,...
                  BasketOptionCall.assetParam.volatility, BasketOptionCall.timeDim.endTime);
              actSolutionCall = genOptPrice(BasketOptionCall);
              BasketOptionCall.payoffParam.putCallType = {'put'};
              BasketOptionPut = BasketOptionCall;
              actSolutionPut = genOptPrice(BasketOptionPut);
              testCase.verifyEqual(actSolutionCall,expSolutionCall);
              testCase.verifyEqual(actSolutionPut,expSolutionPut);
          end
          function testSobolSolution(testCase)
              inp.payoffParam.optType = {'basket'};
              inp.payoffParam.basketWeight = [0.2 0.8];
              inp.assetParam.initPrice = [11 15]; 
              inp.assetParam.volatility = [0.5 0.6];
              inp.assetParam.nAsset = 2;
              inp.assetParam.sqCorr = [1 1; 1 1];
              inp.priceParam.cubMethod = 'Sobol';
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              BasketOptionCall = optPrice(inp);
              [expSolutionCall,expSolutionPut] = vanillaPrice(BasketOptionCall.assetParam.initPrice,...
                  BasketOptionCall.payoffParam.strike, BasketOptionCall.assetParam.interest,...
                  BasketOptionCall.assetParam.volatility, BasketOptionCall.timeDim.endTime);
              actSolutionCall = genOptPrice(BasketOptionCall);
              BasketOptionCall.payoffParam.putCallType = {'put'};
              BasketOptionPut = BasketOptionCall;
              actSolutionPut = genOptPrice(BasketOptionPut);
              testCase.verifyEqual(actSolutionCall,expSolutionCall);
              testCase.verifyEqual(actSolutionPut,expSolutionPut);
          end
      end
      
  end 