inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for three months
inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0.05; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 120; %strike price
inp.payoffParam.putCallType = {'put'}; %looking at a put option
inp.priceParam.absTol = 0.02; %absolute tolerance of a two cents
inp.priceParam.relTol = 0; %zero relative tolerance
EuroPut = optPayoff(inp); %construct an optPrice object 

AmerPut = optPayoff(EuroPut); %construct an American put object
AmerPut.payoffParam.optType = {'american'};

AE = optPayoff(AmerPut);
AE.payoffParam = ...
   struct('putCallType', {{'put','put'}}, ... %note two kinds of option payoffs
   'optType',{{'american','american'}},...
   'strike', [110, 100],...
   'digitalPay', 200) %this needs to have the same dimension
genOptPayoffs(AE,2)


