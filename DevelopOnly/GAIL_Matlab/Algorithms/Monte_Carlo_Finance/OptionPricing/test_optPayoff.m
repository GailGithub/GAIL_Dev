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
   struct('putCallType', {{'put','call'}}, ... %note two kinds of option payoffs
   'optType',{{'digitalcash','digitalcash'}});%,...
   %'strike', 100,...
   %'digitalPay', 200) %this needs to have the same dimension



inp.payoffParam.optType = {'basket'};
inp.payoffParam.basketWeight = [0.2 0.8];
inp.assetParam.nAsset = 2;
inp.assetParam.initPrice = [11 15];
inp.assetParam.corrMat = [1 0.5; 0.5 1];
%Compute the Option Price
BasketOption1 = optPrice(inp);
[Price1, Out1] = genOptPrice(BasketOption1)
     

