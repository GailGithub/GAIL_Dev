%Computing prices of various options
%  Requires OptionOutput.m, payoff.m, stockpath.m,
%  exactprice.m
format compact %remove blank lines from output
clearvars %clear all variables

%% Parameter set-up


inp.timeDim.timeVector=1/64:1/64:0.5; %number of trading periods

inp.assetParam.initPrice=100; %initial stock price
inp.assetParam.volatility=0.5; %volatility
inp.assetParam.interest=0.015; %interest rate


inp.payoffParam.strike=110; %strike price
inp.payoffParam.optType={'american'};
inp.payoffParam.putCallType={'put'};
inp.priceParam.absTol=0.05;
tianci=optPrice(inp)
tic,
[price,out]=genOptPrice(tianci),toc
fred=tianci;
fred.payoffParam.optType={'euro'};
fred.exactPrice