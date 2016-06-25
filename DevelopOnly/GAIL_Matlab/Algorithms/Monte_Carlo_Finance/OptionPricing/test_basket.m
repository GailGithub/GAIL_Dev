inp.payoffParam.optType = {'basket', 'basket','basket','basket'};
inp.payoffParam.putCallType = {'put', 'put','call','call'};
inp.payoffParam.strike = [12 12 14 15];
inp.payoffParam.basketWeight = [0.2 0.8];
inp.assetParam.nAsset = 2;
inp.assetParam.initPrice = [11 15];
inp.assetParam.volatility = [0.5 0.6];
inp.assetParam.corrMat = [1 0.5; 0.5 1];

BO = optPayoff(inp);
genOptPayoffs(BO,2)
