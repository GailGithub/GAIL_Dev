function price=optprice(pay,sample,asset,option)

zvalue=2.58; %for 99% confidence

%% Compute price and error bars for iid sampling w/o control variates
if option.notamer
    price.callN=cumsum(pay.call)./sample.Nvec; %fair price estimated by MC at each N  
end
price.putN=cumsum(pay.put)./sample.Nvec; %fair price estimated by MC at each N
price.callerr='???'; 
price.puterr=price.callerr;

switch sample.type
    case 'iid'
        if option.notamer; 
            if sample.anti==true;
                secall=std((pay.call(1:sample.iidN)+pay.call(sample.iidN+1:sample.N))/2);
            else
                secall=std(pay.call);
            end
            price.callerr=zvalue*secall/sqrt(sample.iidN); 
        end
        if sample.anti==true;
            seput=std((pay.put(1:sample.iidN)+pay.put(sample.iidN+1:sample.N))/2);
        else
            seput=std(pay.put);
        end
        price.puterr=zvalue*seput/sqrt(sample.iidN);
end

%% Compute price and error bars for sampling with control variates
if asset.ncontrol>0 %control variates
    option.exacttype=asset.control;
    controlprice=exactprice(sample,asset,option);
    whegbm=strcmp('eurogbm',asset.control);
    whejump=strcmp('eurojump',asset.control);
    whgmean=strcmp('gmean',asset.control);
    whpriceT=strcmp('priceT',asset.control);
    if any(whegbm)
        callcontrolmean(whegbm)=controlprice.eurogbmcall; 
        putcontrolmean(whegbm)=controlprice.eurogbmput;
    end
    if any(whejump)
        callcontrolmean(whejump)=controlprice.eurojumpcall; 
        putcontrolmean(whejump)=controlprice.eurojumpput;
    end
    if any(whgmean)
        callcontrolmean(whgmean)=controlprice.gmeancall; 
        putcontrolmean(whgmean)=controlprice.gmeanput;   
    end
    if any(whpriceT)
        callcontrolmean(whpriceT)=controlprice.priceT; 
        putcontrolmean(whpriceT)=controlprice.priceT;   
    end
    if option.notamer
        xbarcall=mean(pay.callcontrol,1);
        betacall=(pay.callcontrol-repmat(xbarcall,sample.N,1))\(pay.call-price.callN(sample.N));
        price.callcv=price.callN(sample.N)+(callcontrolmean-xbarcall)*betacall;
    end
    xbarput=mean(pay.putcontrol,1);
    betaput=(pay.putcontrol-repmat(xbarput,sample.N,1))\(pay.put-price.putN(sample.N));
    price.putcv=price.putN(sample.N)+(putcontrolmean-xbarput)*betaput;
    switch sample.type
        case 'iid'
            if option.notamer
                residcall=pay.call-price.callN(sample.N)-(pay.callcontrol-repmat(xbarcall,sample.N,1))*betacall;
                if sample.anti==true;
                    residcall=(residcall(1:sample.iidN)+residcall(sample.iidN+1:sample.N))/2;
                end
                secall=sqrt(sum((residcall).^2,1)/(sample.iidN-asset.ncontrol));
                price.callcverr=zvalue*secall/sqrt(sample.iidN);
            end
            residput=pay.put-price.putN(sample.N)-(pay.putcontrol-repmat(xbarput,sample.N,1))*betaput;
            if sample.anti==true;
                residput=(residput(1:sample.iidN)+residput(sample.iidN+1:sample.N))/2;
            end
            seput=sqrt(sum((residput).^2,1)/(sample.iidN-asset.ncontrol));
            price.putcverr=zvalue*seput/sqrt(sample.iidN);
    end
end
