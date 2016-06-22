inp.payoffParam.optType = {'basket'};
inp.payoffParam.basketWeight = [0.2 0.8];
inp.assetParam.nAsset = 2;
inp.assetParam.initPrice = [11 15];
inp.assetParam.volatility = [0.5 0.6];
inp.assetParam.corrMat = [1 0.5; 0.5 1];
%Compute the Option Price
BasketOption1 = optPrice(inp);
[Price1, Out1] = genOptPrice(BasketOption1)
%Change from call to put
BasketOption2=BasketOption1;
BasketOption2.payoffParam.putCallType = {'put'};
[Price2, Out2] = genOptPrice(BasketOption2)