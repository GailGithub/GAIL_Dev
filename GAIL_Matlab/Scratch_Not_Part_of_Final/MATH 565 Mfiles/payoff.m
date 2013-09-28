function [pay,option]=payoff(smat,sample,asset,option)

%% Payoffs for options
switch option.type
    case 'euro'
        option.name='European';
        Sdminstrike=(smat(:,sample.d+1)-option.strike)*exp(-asset.r*sample.T);
        pay.call=max(Sdminstrike,0); %discounted call payoff
        pay.put=max(-Sdminstrike,0); %discounted put payoff
    case 'amean'
        option.name='Arithmetic Mean Asian';
        Smeanminstrike=(mean(smat(:,2:sample.d+1),2)-option.strike)*exp(-asset.r*sample.T);
        pay.call=max(Smeanminstrike,0); %discounted call payoff
        pay.put=max(-Smeanminstrike,0); %discounted put payoff
    case 'gmean'
        option.name='Geometric Mean Asian';
        Smeanminstrike=(exp(mean(log(smat(:,2:sample.d+1)),2))-option.strike)*exp(-asset.r*sample.T);
        pay.call=max(Smeanminstrike,0); %discounted call payoff
        pay.put=max(-Smeanminstrike,0); %discounted put payoff
    case 'upin'
        option.name='Up and In Barrier';
        Sdminstrike=(smat(:,sample.d+1)-option.strike)*exp(-asset.r*sample.T);
        in=any(smat(:,2:sample.d+1)>option.barrier,2);
        pay.call=max(Sdminstrike.*in,0); %discounted call payoff
        pay.put=max(-Sdminstrike.*in,0); %discounted put payoff
    case 'upout'
        option.name='Up and Out Barrier';
        Sdminstrike=(smat(:,sample.d+1)-option.strike)*exp(-asset.r*sample.T);
        out=any(smat(:,2:sample.d+1)>option.barrier,2);
        pay.call=max(Sdminstrike.*(~out),0); %discounted call payoff
        pay.put=max(-Sdminstrike.*(~out),0); %discounted put payoff
    case 'downin'
        option.name='Down and In Barrier';
        Sdminstrike=(smat(:,sample.d+1)-option.strike)*exp(-asset.r*sample.T);
        in=any(smat(:,2:sample.d+1)<option.barrier,2);
        pay.call=max(Sdminstrike.*(in),0); %discounted call payoff
        pay.put=max(-Sdminstrike.*(in),0); %discounted put payoff
    case 'downout'
        option.name='Down and Out Barrier';
        Sdminstrike=(smat(:,sample.d+1)-option.strike)*exp(-asset.r*sample.T);
        out=any(smat(:,2:sample.d+1)<option.barrier,2);
        pay.call=max(Sdminstrike.*(~out),0); %discounted call payoff
        pay.put=max(-Sdminstrike.*(~out),0); %discounted put payoff
    case 'look'
        option.name='Lookback';
        Sdminstrike=(smat(:,sample.d+1)-min(smat(:,2:sample.d+1),[],2))*exp(-asset.r*sample.T);
        pay.call=max(Sdminstrike,0); %discounted call payoff
        Sdminstrike=(smat(:,sample.d+1)-max(smat(:,2:sample.d+1),[],2))*exp(-asset.r*sample.T);
        pay.put=max(-Sdminstrike,0); %discounted put payoff
    case 'amer'
        option.name='American put';
        basis= @(x) repmat(exp(-x/2),1,3).*[ones(length(x),1) 1-x 1-2*x+x.*x/2];
        sample.tvec=(0:sample.d)*sample.delt; %vector of time values
        putpayoff=max(option.strike-smat,0).*repmat(exp(-asset.r*sample.tvec),sample.N,1); %discounted payoff at each time
        cashflow=putpayoff(:,sample.d+1); %cash flow according to exercise rule
        extime=repmat(sample.d+1,sample.N,1);
        pay.exbound=zeros(sample.d+1,1); pay.exbound(sample.d+1)=option.strike; %exercise boundary
        for j=sample.d:-1:1 %work backwards from expiry to present
            in=find(putpayoff(:,j)>0); %which paths are in the money
            if ~isempty(in)
                if j>1;
                    regmat=basis(smat(in,j)/asset.s0); %regression matrix for stock prices
                    if sample.import==true; %weighted regression for importance sampling
                        [~,numbasis]=size(regmat);
                        regwt=sqrt(sample.impwt(in));
                        hold=regmat*((repmat(regwt,1,numbasis).*regmat)\(regwt.*cashflow(in))); %regressed value of options in the future
                    else
                        hold=regmat*(regmat\cashflow(in)); %regressed value of options in the future
                    end
                else
                    if sample.import==true;
                        hold=mean(cashflow(in).*sample.impwt(in));
                    else
                        hold=mean(cashflow(in));
                    end
                end
                shouldex=in(putpayoff(in,j)>hold); %which paths should be excercised now
                if ~isempty(shouldex);
                    cashflow(shouldex)=putpayoff(shouldex,j); %updated cashflow
                    extime(shouldex)=j;
                    pay.exbound(j)=max(smat(shouldex,j)); 
                end
            end
        end
        pay.put=cashflow;
end

%% Payoffs for control variates
wh=strcmp('eurogbm',asset.control)|strcmp('eurojump',asset.control);
if any(wh);
    Sdminstrike=(smat(:,sample.d+1)-option.strike)*exp(-asset.r*sample.T);
    pay.callcontrol(:,wh)=max(Sdminstrike,0); %discounted call payoff
    pay.putcontrol(:,wh)=max(-Sdminstrike,0); %discounted put payoff
end
wh=strcmp('gmean',asset.control);
if any(wh);
    Smeanminstrike=(exp(mean(log(smat(:,2:sample.d+1)),2))-option.strike)*exp(-asset.r*sample.T);
    pay.callcontrol(:,wh)=max(Smeanminstrike,0); %discounted call payoff
    pay.putcontrol(:,wh)=max(-Smeanminstrike,0); %discounted put payoff
end
wh=strcmp('priceT',asset.control);
if any(wh);
    pay.callcontrol(:,wh)=smat(:,sample.d+1); %price at end time
    pay.putcontrol(:,wh)=smat(:,sample.d+1); %price at end time
end

%% Importance sampling
option.notamer=~strcmp('amer',option.type);
if sample.import==true; 
    if option.notamer; 
        pay.call=pay.call.*sample.impwt;
    end
    pay.put=pay.put.*sample.impwt;
    if asset.ncontrol>0;
        controlwt=repmat(sample.impwt,1,asset.ncontrol);
        if option.notamer; 
            pay.callcontrol=pay.callcontrol.*controlwt;
        end
        pay.putcontrol=pay.putcontrol.*controlwt;
    end
end 

        
