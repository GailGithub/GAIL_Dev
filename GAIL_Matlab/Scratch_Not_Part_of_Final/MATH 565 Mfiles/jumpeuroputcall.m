function [putprice,callprice]=jumpeuroputcall(sample,asset,option)
%normcdf = @(x) 0.5 * erfc(-(x) ./ sqrt(2));
probtol=1e-4;
mp1=exp(asset.ajump+asset.bjump.*asset.bjump/2);
m=mp1-1;
intensity=asset.ljump.*mp1*sample.T;
putprice=0;
callprice=0;
njump=0;
prob=exp(-intensity);
assetn=asset;
dontstop=1;
while dontstop
    assetn.sig=sqrt(asset.sig.*asset.sig+asset.bjump*asset.bjump*njump/sample.T);
    assetn.r=asset.r-asset.ljump*m+njump*log(mp1)/sample.T;
    [gbmput,gbmcall]=europutcall(sample,assetn,option);
    putprice=putprice+prob*gbmput;
    callprice=callprice+prob*gbmcall;
    njump=njump+1;
    oldprob=prob;
    prob=prob*intensity/njump;
    if prob<oldprob && prob<probtol; dontstop=0; end
end
%keyboard
