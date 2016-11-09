%% Importance Sampling
%
% 
%% Set assetPath parameters
clearvars
T=5;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector
inp.assetParam.interest = 0.0;          % risk-free interest rate
inp.assetParam.volatility = 0.3;        % fixed vlatility of asset prices
inp.assetParam.Vinst = 0.09;            % initial value of volatility
inp.assetParam.Vlong = 0.09;            % theta
inp.assetParam.kappa = 1;               % kappa
inp.assetParam.nu = 0;                  % volatility of asset price volatility
inp.assetParam.rho = 0;                 % rho
inp.assetParam.pathType = 'GBM'; 
inp.priceParam.cubMethod = 'IID_MC';
%Set optPayoff parameter
inp.payoffParam.strike = 60;            % strike price
%Set error tolerance
inp.priceParam.absTol = 0;              % absolute tolerance
inp.priceParam.relTol = 0.03;           % three penny on the dollar relative tolerance
initialPrice = [80,60,40];
inp.assetParam.initPrice = 80;         % initial asset price

%% European call option
%With importance sampling
inp.assetParam.meanShift = 5;
ourGBMCallPrice = optPrice(inp)
GBMCallPrice_withIS=zeros(1,3);
out2 = zeros(1,3);
for i=1:3
    inp.assetParam.initPrice = initialPrice(i);
    %Construct an optPrice object
    ourGBMCallPrice = optPrice(inp);
    [GBMCallPrice_withIS(i), out2] = genOptPrice(ourGBMCallPrice);
end
GBMCallPrice_withIS
%Without importance Sampling
inp.assetParam.meanShift = 0;
GBMCallPrice_withoutIS=zeros(1,3);
out1 = zeros(1,3);
for i=1:3
    inp.assetParam.initPrice = initialPrice(i);
    %Construct an optPrice object
    ourGBMCallPrice = optPrice(inp);
    [GBMCallPrice_withoutIS(i), out1] = genOptPrice(ourGBMCallPrice);
end
GBMCallPrice_withoutIS

%% European put option
inp.payoffParam.putCallType = {'put'};
initialPrice = fliplr(initialPrice);   %initialPrice = [40,60,80]
% With Importance Sampling
inp.assetParam.meanShift = -2;
GBMPutPrice_withIS=zeros(1,3);
out4 = zeros(1,3);
ourGBMPutPrice = optPrice(inp)
for i=1:3
    inp.assetParam.initPrice = initialPrice(i);
    %Construct an optPrice object
    ourGBMPutPrice = optPrice(inp);
    [GBMPutPrice_withIS(i), out4] = genOptPrice(ourGBMPutPrice);
end
GBMPutPrice_withIS
%Without Importance Sampling
inp.assetParam.meanShift = 0;
GBMPutPrice_withoutIS=zeros(1,3);
out3 = zeros(1,3);
for i=1:3
    inp.assetParam.initPrice = initialPrice(i);
    %Construct an optPrice object
    ourGBMPutPrice = optPrice(inp);
    [GBMPutPrice_withoutIS(i), out3] = genOptPrice(ourGBMPutPrice);
end
GBMPutPrice_withoutIS

%% American put option
inp.payoffParam.optType = {'american'}; %change from European to American
inp.assetParam.meanShift = -2;
AmericanPutPrice_withIS=zeros(1,3);
out6 = zeros(1,3);
AmericanPut = optPrice(inp)
for i=1:3
    inp.assetParam.initPrice = initialPrice(i);
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    [AmericanPutPrice_withIS(i), out6] = genOptPrice(AmericanPut);
end
AmericanPutPrice_withIS

inp.assetParam.meanShift = 0;
AmericanPutPrice_withoutIS=zeros(1,3);
out5 = zeros(1,3);
for i=1:3
    inp.assetParam.initPrice = initialPrice(i);
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    [AmericanPutPrice_withoutIS(i), out5] = genOptPrice(AmericanPut);
end
AmericanPutPrice_withoutIS
%% Reference

