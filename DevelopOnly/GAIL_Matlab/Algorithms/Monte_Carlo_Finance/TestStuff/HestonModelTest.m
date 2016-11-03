classdef HestonModelTest < matlab.unittest.TestCase
      % Heston Model Testing 
      % Compared with geometric Brownian Motion.
      % Authors: Tianci Zhu
    
      methods (Test)
          % Test Heston Model using Black Scholes formula
          function testIIDSolution(testCase)
              % Assign parameter values for european option
              euroinp.payoffParam.optType = {'euro'};
              T=1;
              delta_t=0.1;
              t0 = delta_t;
              euroinp.timeDim.timeVector = t0:delta_t:T;
              euroinp.assetParam.interest=0;
              euroinp.assetParam.kappa = 1;
              euroinp.assetParam.nu = 0;%1e-16;
              euroinp.priceParam.absTol = 0;
              euroinp.priceParam.relTol = 0.01;
              euroinitPrice=[100 120];
              eurovolatility=[0.5 0.3];
              strike=[90 110];
              rho=[0.3 0.8];
              callPut=[{'call'},{'put'}];
              nOption=2;
              q=1;
              QEPrice=zeros(1,2*nOption^5);
              QEmPrice=zeros(1,2*nOption^5);
              exactPrice=zeros(1,2*nOption^5);
              for i=1:nOption
                  for j=1:nOption
                      for m=1:nOption
                          for n=1:2
                              for v=1:nOption
                                  for r=1:nOption
                                      %euroinp.assetParam.pathType='GBM';
                                      euroinp.assetParam.pathType = 'QE';
                                      euroinp.payoffParam.putCallType = callPut(n);
                                      euroinp.payoffParam.strike = strike(m);
                                      euroinp.assetParam.volatility = eurovolatility(j);
                                      euroinp.assetParam.initPrice = euroinitPrice(i);
                                      euroinp.assetParam.Vinst = eurovolatility(j)^2; 
                                      euroinp.assetParam.Vlong = eurovolatility(j)^2;
                                      euroinp.assetParam.rho = rho(r);
                                      eurooption = optPrice(euroinp);
                                      QEPrice(q)=genOptPrice(eurooption);
                                      euroinp.assetParam.pathType = 'QE_m';
                                      eurooption = optPrice(euroinp);
                                      QEmPrice(q)=genOptPrice(eurooption);
                                      exactPrice(q)=eurooption.exactPrice;
                                      q = q+1;
                                  end
                              end
                          end
                      end
                  end
              end
              testCase.verifyLessThan(abs(QEPrice-exactPrice)./exactPrice,eurooption.priceParam.relTol);
              testCase.verifyLessThan(abs(QEmPrice-exactPrice)./exactPrice,eurooption.priceParam.relTol);
              
          end            
          
     end
      
end 