  classdef BasketOptionTest < matlab.unittest.TestCase
      % BasketOptionTest tests the basket option price generated using
      % Monte Carlo Method and Quasi Monte Carlo Method.
    
      methods (Test)
          function testIIDSolution(testCase)
              euroinp.timeDim.timeVector = 1;
              euroinp.payoffParam.optType = {'euro'};
              euroinp.payoffParam.putCallType = {'call'};
              eurooption = optPayoff(euroinp);
              exactcallprice=eurooption.exactPrice;
              eurooption.payoffParam.putCallType = {'put'};
              exactputprice=eurooption.exactPrice;
              inp.payoffParam.optType = {'basket'};
              inp.payoffParam.basketWeight = [0.3 0.7];
              inp.assetParam.initPrice = [10 10]; 
              inp.assetParam.volatility = [0.5 0.5];
              inp.assetParam.nAsset = 2;
              inp.assetParam.corrMat = [1 1;1 1];
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              BasketOption = optPrice(inp);
               
               appSolutionCall = genOptPrice(BasketOption);
               BasketOption.payoffParam.putCallType = {'put'};
               appSolutionPut = genOptPrice(BasketOption);
               testCase.verifyLessThan(abs(appSolutionCall-exactcallprice)/exactcallprice,BasketOption.priceParam.relTol);
               testCase.verifyLessThan(abs(appSolutionPut-exactputprice)/exactputprice,BasketOption.priceParam.relTol);
           end

          function testSobolSolution(testCase)
              euroinp.timeDim.timeVector = 1;
              euroinp.payoffParam.optType = {'euro'};
              euroinp.payoffParam.putCallType = {'call'};
              euroinp.priceParam.cubMethod = 'Sobol';
              eurooption = optPayoff(euroinp);
              exactcallprice=eurooption.exactPrice;
              eurooption.payoffParam.putCallType = {'put'};
              exactputprice=eurooption.exactPrice;
              basketinp.payoffParam.optType = {'basket'};
              basketinp.payoffParam.basketWeight = [0.3 0.7];
              basketinp.assetParam.initPrice = [10 10]; 
              basketinp.assetParam.volatility = [0.5 0.5];
              basketinp.assetParam.nAsset = 2;
              basketinp.assetParam.corrMat = [1 1;1 1];
              basketinp.priceParam.absTol = 0;   
              basketinp.priceParam.relTol = 0.01;
              basketinp.priceParam.cubMethod = 'Sobol';
              basketinp.timeDim.timeVector = 1;
              
              BasketOption = optPrice(basketinp);
              appSolutionCall = genOptPrice(BasketOption);
              BasketOption.payoffParam.putCallType = {'put'};
              appSolutionPut = genOptPrice(BasketOption);
              testCase.verifyLessThan(abs(appSolutionCall-exactcallprice)/exactcallprice,BasketOption.priceParam.relTol);
              testCase.verifyLessThan(abs(appSolutionPut-exactputprice)/exactputprice,BasketOption.priceParam.relTol);
          end
      end
      
  end 