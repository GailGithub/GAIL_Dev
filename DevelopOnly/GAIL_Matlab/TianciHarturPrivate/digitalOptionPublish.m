%% Digital Option in Matlab
% * An option whose payout is fixed after the underlying stock exceeds the predetermined threshold or strike price. 


%Pricing digital option with cash
inp.payoffParam.optType={'digital'};
inp.priceParam.absTol=0.01;
tianci=optPrice(inp)
hartur=genOptPrice(tianci)

%Pricing digital option with asset
inp.payoffParam.optType={'digital'};
inp.payoffParam.cashAssetType={'asset'};
inp.priceParam.absTol=0.01;
tianci=optPrice(inp)
hartur=genOptPrice(tianci)
