function price=exactprice(sample,asset,option)

%Expected value of ending price of asset
if any(strcmp('priceT',option.exacttype)); 
    price.priceT=asset.s0*exp(asset.r*sample.T);
end

%Pricing European geometric brownian motion
if any(strcmp('eurogbm',option.exacttype)); 
    [price.eurogbmcall,price.eurogbmput]=eurogbmprice(sample,asset,option);
end

%Pricing European geometric brownian motion with jumps
if any(strcmp('eurojump',option.exacttype)); 
    probtol=1e-4;
    mp1=exp(asset.ajump+asset.bjump.*asset.bjump/2);
    m=mp1-1;
    intensity=asset.ljump.*mp1*sample.T;
    price.eurojumpcall=0;
    price.eurojumpput=0;
    njump=0;
    prob=exp(-intensity);
    assetn=asset;
    dontstop=1;
    while dontstop
        assetn.sig=sqrt(asset.sig.*asset.sig+asset.bjump*asset.bjump*njump/sample.T);
        assetn.r=asset.r-asset.ljump*m+njump*log(mp1)/sample.T;
        [gbmcall,gbmput]=eurogbmprice(sample,assetn,option);
        price.eurojumpcall=price.eurojumpcall+prob*gbmcall;
        price.eurojumpput=price.eurojumpput+prob*gbmput;
        njump=njump+1;
        oldprob=prob;
        prob=prob*intensity/njump;
    if prob<oldprob && prob<probtol; dontstop=0; end
    end
end

%Pricing Asian geometric mean
if any(strcmp('gmean',option.exacttype)); 
    assetbar=asset; samplebar=sample;
    samplebar.T=(1+1/sample.d)*sample.T/2; 
    assetbar.sig=asset.sig*sqrt((2+1/sample.d)/3);
    assetbar.r=asset.r+(assetbar.sig^2-asset.sig^2)/2;
    [price.gmeancall,price.gmeanput]=eurogbmprice(samplebar,assetbar,option);
    price.gmeancall=price.gmeancall*exp(assetbar.r*samplebar.T - asset.r*sample.T);
    price.gmeanput=price.gmeanput*exp(assetbar.r*samplebar.T - asset.r*sample.T);
end


function [callprice,putprice]=eurogbmprice(sample,asset,option)
%normcdf = @(x) 0.5 * erfc(-(x) ./ sqrt(2));
priceratio=option.strike*exp(-asset.r*sample.T)./asset.s0;
xbig=log(priceratio)./(asset.sig*sqrt(sample.T))+asset.sig*sqrt(sample.T)/2;
xsmall=log(priceratio)./(asset.sig*sqrt(sample.T))-asset.sig*sqrt(sample.T)/2;
putprice=asset.s0.*(priceratio.*normcdf(xbig) - normcdf(xsmall));
callprice=putprice+asset.s0*(1-priceratio);