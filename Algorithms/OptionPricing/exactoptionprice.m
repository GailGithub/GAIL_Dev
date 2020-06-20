function price=exactoptionprice(stparam,payparam,optiontype)
price=zeros(1,numel(optiontype));

%Expected value of ending price of asset
wh=strcmp('stockprice',optiontype);
if any(wh)
   price(wh)=stparam.S0;
end

%Pricing European geometric brownian motion
whcall=strcmp('eurocall',optiontype);
whput=strcmp('europut',optiontype);
if any(whcall | whput) 
   [eurocall,europut]=eurogbmprice(stparam,payparam);
   price(whcall)=eurocall;
   price(whput)=europut;
end

%Pricing Asian geometric mean
whcall=strcmp('gmeancall',optiontype);
whput=strcmp('gmeanput',optiontype);
if any(whcall | whput)
   stparambar=stparam;
   stparambar.T=(1+1/stparam.d)*stparam.T/2; 
   stparambar.sig=stparam.sig*sqrt((2+1/stparam.d)/3);
   stparambar.r=stparam.r+(stparambar.sig^2-stparam.sig^2)/2;
   [gmeancall,gmeanput]=eurogbmprice(stparambar,payparam);
   gmeancall=gmeancall*exp(stparambar.r*stparambar.T - stparam.r*stparam.T);
   gmeanput=gmeanput*exp(stparambar.r*stparambar.T - stparam.r*stparam.T);
   price(whcall)=gmeancall;
   price(whput)=gmeanput;
end

function [callprice,putprice]=eurogbmprice(stparam,payparam)
%normcdf = @(x) 0.5 * erfc(-(x) ./ sqrt(2));
priceratio=payparam.K*exp(-stparam.r*stparam.T)./stparam.S0;
xbig=log(priceratio)./(stparam.sig*sqrt(stparam.T))+stparam.sig*sqrt(stparam.T)/2;
xsmall=log(priceratio)./(stparam.sig*sqrt(stparam.T))-stparam.sig*sqrt(stparam.T)/2;
putprice=stparam.S0.*(priceratio.*normcdf(xbig) - normcdf(xsmall));
callprice=putprice+stparam.S0*(1-priceratio);