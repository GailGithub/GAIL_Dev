 classdef americanOptionTest < matlab.unittest.TestCase
      % americanOptionTest tests the American option price generated using
      % Longstaff Least Quare Method
      % Authors: Tianci Zhu
    
    methods (Test)
          % Test Monte Carlo Method
          function testIIDSolution(testCase)
              % Computing prices of various options
              % Requires OptionOutput.m, payoff.m, stockpath.m,
              % exactprice.m
              format compact %remove blank lines from output
              clearvars %clear all variables

              %% Parameter set-up
              inp.payoffParam.optType={'american'};
              inp.payoffParam.putCallType={'put'};
              inp.assetParam.interest=0.06; %interest rate
              inp.payoffParam.strike=40; %strike price
              inp.priceParam.absTol=0.05;
              initPrice=[36 38 40 42 44];
              year=[1 2];
              volatility=[0.2 0.4];
              p=1;
              for i=1:5
                  for j=1:2
                      for m=1:2
                      inp.assetParam.initPrice =initPrice(i); %initial stock price
                      inp.timeDim.timeVector=1/50:1/50:year(j); %number of trading periods
                      inp.assetParam.volatility=volatility(m); %volatility
                      americanOption=optPrice(inp);
                      actualPrice(p)=genOptPrice(americanOption);
                      p=p+1;
                      end
                  end
              end
            
              %% Enter the value given in Longstaff's paper
              longstaffValue=[4.472 7.091 4.821 8.488 3.244 6.139 3.735... 
                  7.669 2.313 5.308 2.879 6.921 1.617 4.588 2.206 6.243...
                  1.118 3.957 1.675 5.622];
              
              %% Test the absolute error
              testCase.verifyLessThan(abs(longstaffValue-actualPrice),americanOption.priceParam.absTol);
          end
    end
 end
 
          

