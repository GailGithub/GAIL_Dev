%Option Price Output
close all
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18)
notamer=~strcmp('amer',option.type);

%% Plot stock paths
if output.stockpath
    numplot=200;
    tplot=(0:sample.d)*sample.delt*365;
    h=plot(tplot,smat(1:min(sample.N,numplot),:),'-',[tplot(1) tplot(sample.d+1)],[option.strike option.strike],'r--'); set(h,'linewidth',2)
    xlabel('Time (days)')
    ylabel('Stock Price')
    eval(['print -depsc MCStockN' int2str(sample.N) 'sample.d' int2str(sample.d) '.eps'])
    refresh
end

%% Plot MC method estimates converging to the answer
if output.MCconverge
    if notamer;
        numplot=200;
        Nplot=floor(10.^((0:numplot)'*(log10(sample.N)/numplot)));
        figure; 
        h=semilogx(Nplot,price.callN(Nplot),'-b',Nplot,price.putN(Nplot),'r-');  
        set(h,'linewidth',2);
        xlabel('Number of Stock Paths')
        ylabel('Option Price')
        legend('Call','Put','Location','Best')
        set(gca,'Xtick',10.^(1:floor(log10(sample.N))))
        eval(['print -depsc Conv' option.type 'Price' int2str(sample.N) '.eps'])
        refresh
    else
        okexbound=pay.exbound>0;
        figure; 
        h=plot(sample.tvec(okexbound),pay.exbound(okexbound),'b-',[0 sample.T],[option.strike option.strike],'r--');  
        set(h,'linewidth',2);
        axis([0 sample.T 0.9*min(pay.exbound(okexbound)) option.strike*1.1])
        xlabel('Time')
        ylabel('Option Price')
        legend('Exercise Boundary','Strike','Location','Best')
        eval(['print -depsc ExBoundAmerPut' int2str(sample.N) '.eps'])
        refresh
    end
end

%% Output numerical results
if output.numer
    disp(['Using ' int2str(sample.N) ' asset price samples based on a'])
    disp(['   ' asset.modname ' model of the asset'])
    disp('   with sampling method:')
    for i=1:numel(sample.name); disp(sample.name{i}); end
    disp(['For an initial asset price of $' num2str(asset.s0,'%0.2f')])
    disp(['   a strike price of $' num2str(option.strike,'%0.2f')])
    if any(strcmp(option.type,{'upin','upout','downin','downout'}))
        disp(['   a barrier of $' num2str(option.barrier,'%0.2f')])
    end
    disp(['   ' num2str(sample.T,'%0.2f') ' years to maturity'])
    disp(['   an interest rate of ' num2str(100*asset.r,'%0.2f') '%'])
    disp(['   a volatility of ' num2str(100*asset.sig,'%0.2f') '%:'])
    disp(['For ' option.name ' options monitored ' int2str(sample.d) ' times'])
    switch asset.modtype 
        case 'stickydelta'
            disp(['   where the skew coefficient is ' num2str(asset.sigskew,'%0.2f') ])
            disp(['   and the smile coefficient is ' num2str(asset.sigsmile,'%0.2f') ])
        case 'jump'
            disp(['   where lognormal jumps occur on average ' num2str(asset.ljump,'%0.2f')  ' times per year'])
            disp(['   and these jumps have mean ' num2str(asset.ajump,'%0.2f') ' and standard deviation ' num2str(asset.bjump,'%0.2f')])
        case 'vg'
            disp(['   with parameter beta = ' num2str(asset.vgbeta,'%0.4f')])
    end
    if asset.ncontrol==0;
        if notamer; 
            disp(['The call price is $' num2str(price.callN(sample.N),'%0.4f') char(177) num2str(price.callerr,'%0.4f')])
            disp(['   and the put price is $' num2str(price.putN(sample.N),'%0.4f') char(177) num2str(price.puterr,'%0.4f')])
        else
            disp(['The put price is $' num2str(price.putN(sample.N),'%0.4f') char(177) num2str(price.puterr,'%0.4f')])
        end
    else
        if notamer; 
            disp('Without control variates:')
            disp(['   the call price is $' num2str(price.callN(sample.N),'%0.4f') char(177) num2str(price.callerr,'%0.4f')])
            disp(['   and the put price is $' num2str(price.putN(sample.N),'%0.4f') char(177) num2str(price.puterr,'%0.4f')])
            disp(['With control variates ' asset.control{:}])
            disp(['   the call price is $' num2str(price.callcv,'%0.4f') char(177) num2str(price.callcverr,'%0.4f')])
            disp(['   and the put price is $' num2str(price.putcv,'%0.4f') char(177) num2str(price.putcverr,'%0.4f')])
        else
            disp('Without control variates:')
            disp(['   the put price is $' num2str(price.putN(sample.N),'%0.4f') char(177) num2str(price.puterr,'%0.4f')])
            disp(['With control variates ' asset.control{:}])
            disp(['   the put price is $' num2str(price.putcv,'%0.4f') char(177) num2str(price.putcverr,'%0.4f')])
        end
    end

    if any(strcmp('eurogbm',option.exacttype));
        disp(['Compared to the GBM European call price of $' num2str(exprice.eurogbmcall,'%0.4f')])
        disp(['   and the GBM European put price of $' num2str(exprice.eurogbmput,'%0.4f')])
    end
    if any(strcmp('eurojump',option.exacttype));
        disp(['Compared to the jump diffusion European call price of $' num2str(exprice.eurojumpcall,'%0.4f')])
        disp(['   and the jump diffusion European put price of $' num2str(exprice.eurojumpput,'%0.4f')])
    end
    if any(strcmp('gmean',option.exacttype));
        disp(['Compared to the Asian geometric mean call price of $' num2str(exprice.gmeancall,'%0.4f')])
        disp(['   and the Asian geometric mean put price of $' num2str(exprice.gmeanput,'%0.4f')])
    end
    disp(['This computation took ' num2str(timetaken) ' seconds'])
    disp(' ')
end
