%% Perform option pricing with meanMC_CLT
%This MATLAB script shows how to use approximate Central Limit Theorem
% (CLT) with Monte Carlo and Control Variate to price a financial
% derivative or option. 

%% Initialize option parameters for a European call option
inp.timeDim.timeVector = 0.25:0.25:1; %time increments of 0.004 up to time 1
inp.assetParam.initPrice = 120; %initial stock price
inp.assetParam.interest = 0.01; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 130; %strike price
inp.priceParam.absTol = 0.1; %absolute tolerance
inp.priceParam.relTol = 0; %relative tolerance
inp.priceParam.cubMethod = 'IID_MC_CLT'; %use MC CLT
EuroCall = optPrice(inp);


%% Price the European Call
[priceEuroCall, outEuro] = genOptPrice(EuroCall)

%% Create an Asian Arithmetic Mean Call option and Price it
ArithMeanCall = optPrice(EuroCall);
ArithMeanCall.payoffParam.optType = {'amean'}
[priceAmeanCall, outAmean] = genOptPrice(ArithMeanCall)


%% Create an Asian Arithmetic Mean Call option with Euro as control variate and Price it
ArithMeanEuroCVCall = optPrice(EuroCall);
ArithMeanEuroCVCall.payoffParam = ...
    struct('optType', {{'amean','gmean','euro'}}, 'putCallType', {{'call','call','call'}})
[priceAmeanEuroCVCall, outEuroCVAmean] = genOptPrice(ArithMeanEuroCVCall)

return
 
 
   
            
   

            
            
