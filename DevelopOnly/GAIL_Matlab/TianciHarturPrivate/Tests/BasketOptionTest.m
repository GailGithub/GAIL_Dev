  classdef BasketOptionTest < matlab.unittest.TestCase
      % BasketOptionTest tests the basket option price generated using
      % Monte Carlo Method and Quasi Monte Carlo Method.
    
      methods (Test)
          function testIIDSolution(testCase)
              euroinp.timeDim.timeVector = 1;
              euroinp.payoffParam.optType = {'euro'};
              eurooption = optPayoff(euroinp);
              exactcallprice=eurooption.exactPrice;
              inp.payoffParam.optType = {'basket'};
              inp.payoffParam.basketWeight = [0.5 0.5];
              inp.assetParam.initPrice = [10 10]; 
              inp.assetParam.volatility = [0.5 0.5];
              inp.assetParam.nAsset = 2;
              inp.assetParam.corrMat = [1 1;1 1];
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              BasketOption = optPrice(inp);
               
               actSolutionCall = genOptPrice(BasketOption);
               BasketOption.payoffParam.putCallType = {'put'};
               %actSolutionPut = genOptPrice(BasketOption);
               testCase.verifyLessThan(abs(actSolutionCall-exactcallprice)/exactcallprice,BasketOption.priceParam.relTol);
               %testCase.verifyLessThan(abs(actSolutionPut-expSolutionPut)/expSolutionPut,0.01);
           end

          function testSobolSolution(testCase)
              euroinp.timeDim.timeVector = 1;
              euroinp.payoffParam.optType = {'euro'};
              euroinp.priceParam.cubMethod = 'Sobol';
              eurooption = optPayoff(euroinp);
              exactcallprice=eurooption.exactPrice;
              basketinp.payoffParam.optType = {'basket'};
              basketinp.payoffParam.basketWeight = [0.5 0.5];
              basketinp.assetParam.initPrice = [10 10]; 
              basketinp.assetParam.volatility = [0.5 0.5];
              basketinp.assetParam.nAsset = 2;
              basketinp.assetParam.corrMat = [1 1;1 1];
              basketinp.priceParam.absTol = 0;   
              basketinp.priceParam.relTol = 0.01;
              basketinp.priceParam.cubMethod = 'Sobol';
              basketinp.timeDim.timeVector = 1;
              
              BasketOption = optPrice(basketinp);
              actSolutionCall = genOptPrice(BasketOption);
              BasketOption.payoffParam.putCallType = {'put'};
              %actSolutionPut = genOptPrice(BasketOption);
              testCase.verifyLessThan(abs(actSolutionCall-exactcallprice)/exactcallprice,BasketOption.priceParam.relTol);
              %testCase.verifyLessThan(abs(actSolutionPut-expSolutionPut)/expSolutionPut,basketOption.priceParam.relTol);
          end
      end
      
  end 