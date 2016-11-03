  classdef BasketOptionTest < matlab.unittest.TestCase
      % BasketOptionTest tests the basket option price generated using
      % Monte Carlo Method and Quasi Monte Carlo Method.
      % Authors: Tianci Zhu and Hartur Santi
    
      methods (Test)
          % Test Monte Carlo Method
          function testIIDSolution(testCase)
              % Assign parameter values for european option
              euroinp.payoffParam.optType = {'euro'};
              euroinp.timeDim.timeVector = 1;
              euroinp.priceParam.absTol = 0;
              euroinp.priceParam.relTol = 0.01;
              euroinp.timeDim.timeVector = 1;
              euroinitPrice=[10 11 12];
              eurovolatility=[0.5 0.6 0.3];
              strike=[9 10 11];
              callPut=[{'call'},{'put'}];
              nOption=length(euroinitPrice);
              q=0;
              for i=1:nOption
                  for j=1:nOption
                      for m=1:nOption
                          for n=1:2
                              for v=1:nOption
                              euroinp.payoffParam.putCallType = callPut(n);
                              euroinp.payoffParam.strike = strike(m);
                              euroinp.assetParam.volatility = eurovolatility(j);
                              euroinp.assetParam.initPrice = euroinitPrice(i);
                              eurooption = optPayoff(euroinp);
                              q = q+1;
                              exactPrice(q)=eurooption.exactPrice;
                              end
                          end
                      end
                  end
              end
              % Set parameter value for basket option
              inp.payoffParam.optType = {'basket'};
              inp.assetParam.nAsset = 2;
              inp.assetParam.corrMat = [1 1;1 1];
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              basketinitPrice = repmat(euroinitPrice',1,inp.assetParam.nAsset);
              basketvolatility = repmat(eurovolatility',1,inp.assetParam.nAsset);
              basketWeight = [0.3 0.7; 0.5 0.5; 0.2 0.8];
              p=0;
              for a=1:nOption
                  for b=1:nOption
                      for c=1:nOption
                          for d=1:2
                              for e=1:nOption
                                 inp.payoffParam.basketWeight = basketWeight(e,:);
                                 inp.payoffParam.putCallType = callPut(d);
                                 inp.payoffParam.strike = strike(c);
                                 inp.assetParam.volatility = basketvolatility(b,:);
                                 inp.assetParam.initPrice = basketinitPrice(a,:);
                                 p=p+1;
                                 BasketOption = optPrice(inp);
                                 % calculate basket option
                                 appPrice(p) = genOptPrice(BasketOption);
                              end
                          end
                      end
                  end
              end
              
              % Test the relative error
              testCase.verifyLessThan(abs(appPrice-exactPrice)./exactPrice,BasketOption.priceParam.relTol);
          end
          
            % Test Quasi-Monte Carlo Method
          function testSobolSolution(testCase)
              % Assign parameter values for european option
              euroinp.payoffParam.optType = {'euro'};
              euroinp.timeDim.timeVector = 1;
              euroinp.priceParam.absTol = 0;
              euroinp.priceParam.relTol = 0.01;
              euroinp.timeDim.timeVector = 1;
              euroinitPrice=[10 11 12];
              eurovolatility=[0.5 0.6 0.3];
              strike=[9 10 11];
              callPut=[{'call'},{'put'}];
              nOption=length(euroinitPrice);
              q=0;
              for i=1:nOption
                  for j=1:nOption
                      for m=1:nOption
                          for n=1:2
                              for v=1:nOption
                              euroinp.payoffParam.putCallType = callPut(n);
                              euroinp.payoffParam.strike = strike(m);
                              euroinp.assetParam.volatility = eurovolatility(j);
                              euroinp.assetParam.initPrice = euroinitPrice(i);
                              eurooption = optPayoff(euroinp);
                              q = q+1;
                              exactPrice(q)=eurooption.exactPrice;
                              end
                          end
                      end
                  end
              end
              % Set parameter value for basket option
              inp.payoffParam.optType = {'basket'};
              inp.assetParam.nAsset = 2;
              inp.assetParam.corrMat = [1 1;1 1];
              inp.priceParam.cubMethod='Sobol';
              inp.priceParam.absTol = 0;   
              inp.priceParam.relTol = 0.01;
              inp.timeDim.timeVector = 1;
              basketinitPrice = repmat(euroinitPrice',1,inp.assetParam.nAsset);
              basketvolatility = repmat(eurovolatility',1,inp.assetParam.nAsset);
              basketWeight = [0.3 0.7; 0.5 0.5; 0.2 0.8];
              p=0;
              for a=1:nOption
                  for b=1:nOption
                      for c=1:nOption
                          for d=1:2
                              for e=1:nOption
                                 inp.payoffParam.basketWeight = basketWeight(e,:);
                                 inp.payoffParam.putCallType = callPut(d);
                                 inp.payoffParam.strike = strike(c);
                                 inp.assetParam.volatility = basketvolatility(b,:);
                                 inp.assetParam.initPrice = basketinitPrice(a,:);
                                 p=p+1;
                                 BasketOption = optPrice(inp);
                                 % calculate basket option
                                 appPrice(p) = genOptPrice(BasketOption);
                              end
                          end
                      end
                  end
              end
              % Test the relative error
              testCase.verifyLessThan(abs(appPrice-exactPrice)./exactPrice,BasketOption.priceParam.relTol);
          end
          
     end
      
  end 