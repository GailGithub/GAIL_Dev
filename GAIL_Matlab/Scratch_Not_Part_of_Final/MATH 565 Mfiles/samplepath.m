function [smat,sample,asset]=samplepath(sample,asset,option)

if sample.anti==true; sample.iidN=sample.N/2; else sample.iidN=sample.N; end

switch sample.type
    case 'iid'
        sample.name={'     independent and identically distributed'};
        if sample.anti==true; sample.name=[sample.name '     with antithetic variates']; end
        if any(strcmp(asset.modtype,{'dgbm','stickydelta','jump','vg' }))
            x=zeros(sample.N,sample.d);
            x(1:sample.iidN,:)=randn(sample.iidN,sample.d); %normal random numbers
            if sample.anti==true; x=[x(1:sample.iidN,:);-x(1:sample.iidN,:)]; end %antithetic sampling
        end
        if any(strcmp(asset.modtype,{'jump','vg' }))
            z=zeros(sample.N,sample.d);
            z(1:sample.iidN,:)=rand(sample.iidN,sample.d); %uniform random numbers
            if sample.anti==true; z=[z(1:sample.iidN,:);1-z(1:sample.iidN,:)]; end %antithetic sampling
        end
        switch asset.modtype
            case 'jump'
                y=zeros(sample.N,sample.d);
                y(1:sample.iidN,:)=randn(sample.iidN,sample.d); %normal random numbers
                if sample.anti==true; y=[y(1:sample.iidN);-y(1:sample.iidN)]; end %antithetic sampling
                z=poissinv(z,asset.ljump*sample.delt); %Poisson random numbers, number of jumps
            case 'vg'
                z=gaminv(z,sample.delt/asset.vgbeta); %Gamma random numbers, size of "Delta"
            case 'cgbm'
                x=zeros(sample.N,sample.s);
                x(1:sample.iidN,:)=randn(sample.iidN,sample.s); %normal random numbers
                if sample.anti==true; x=[x(1:sample.iidN,:);-x(1:sample.iidN,:)]; end %antithetic sampling
        end
    case 'sobol'
        sample.name={'     Sobol quasirandom'};
        if sample.anti==true; sample.name=[sample.name '     with antithetic variates']; end
        if any(strcmp(asset.modtype,{'dgbm'}))
            pscsob=scramble(sobolset(sample.d),'MatousekAffineOwen');
            x=zeros(sample.N,sample.d);
            x(1:sample.iidN,:)=norminv(net(pscsob,sample.iidN)); %normal random numbers
            if sample.anti==true; x=[x(1:sample.iidN,:);-x(1:sample.iidN,:)]; end %antithetic sampling
        end
        switch asset.modtype
            case {'KL','BB'}
                pscsob=scramble(sobolset(sample.s),'MatousekAffineOwen');
                x=zeros(sample.N,sample.s);
                x(1:sample.iidN,:)=norminv(net(pscsob,sample.iidN)); %normal random numbers
                if sample.anti==true; x=[x(1:sample.iidN,:);-x(1:sample.iidN,:)]; end %antithetic sampling
        end       
end

if sample.import==true %importance sampling, mean shift
    sample.name=[sample.name '     with importance sampling'];
    sample.name=[sample.name ['         mean shift = ' num2str(sample.meanshift)]];
    xd=size(x,2);
    if numel(sample.meanshift)~=xd; sample.meanshift=repmat(sample.meanshift(1),1,xd); end
    x=x+repmat(sample.meanshift,sample.N,1);
    sample.impwt=exp(sample.meanshift*sample.meanshift'/2-x*sample.meanshift');
end


asset.ncontrol=numel(asset.control);
switch asset.modtype
    case 'dgbm'
        asset.modname='discrete geometric Brownian motion';
        smat=cumprod([asset.s0*ones(sample.N,1) exp((asset.r-asset.sig.*asset.sig/2)*sample.delt + x*asset.sig*sqrt(sample.delt))],2);
    case 'stickydelta'
        asset.modname='DGBM with deterministic volatility skew and smile';
        smat=[asset.s0*ones(sample.N,1) zeros(sample.N,sample.d)];
        for j=1:sample.d
            vol=asset.sig+asset.sigskew*(smat(:,j)./option.strike-1)+asset.sigsmile*(smat(:,j)./option.strike-1).^2;
            smat(:,j+1)=smat(:,j).*exp((asset.r-vol.*vol/2)*sample.delt + x(:,j).*vol*sqrt(sample.delt));
        end
    case 'jump'
        asset.modname='DGBM with jumps';
        smat=zeros(sample.N,sample.d+1);
        mu=asset.r-asset.sig*asset.sig/2-asset.ljump.*(exp(asset.ajump+asset.bjump.*asset.bjump/2)-1);
        for j=1:sample.d
            smat(:,j+1)=smat(:,j)+mu*sample.delt + x(:,j).*(asset.sig*sqrt(sample.delt))+asset.ajump*z(:,j)+asset.bjump*sqrt(z(:,j)).*y(:,j);
        end
        smat=asset.s0*exp(smat);
    case 'vg'
        asset.modname='variance gamma';
        smat=zeros(sample.N,sample.d+1);
        mu=asset.r+log(1-asset.vgbeta*asset.sig*asset.sig/2)/asset.vgbeta;
        for j=1:sample.d
            smat(:,j+1)=smat(:,j)+mu*sample.delt + (asset.sig*sqrt(asset.vgbeta)).*x(:,j).*sqrt(z(:,j));
        end
        smat=asset.s0*exp(smat);
    case 'KL'
        asset.modname=['continuous geometric Brownian motion (Karhunen-Loeve) with ' int2str(sample.s) ' eigenfunctions'];
        tovT=(0:sample.d)/sample.d;
        xmat=zeros(sample.N,sample.d+1);
        km1=pi*((1:sample.s)-1/2);
        for k=1:sample.s;
            xmat=xmat+x(:,k)*(sin(km1(k)*tovT)*sqrt(2*sample.T)/km1(k));
        end
        smat=asset.s0*exp((asset.r-asset.sig.*asset.sig/2)*repmat(tovT*sample.T,sample.N,1) + xmat*asset.sig);
    case 'BB'
        asset.modname=['continuous geometric Brownian motion (Brownian bridge) with ' int2str(sample.s) ' eigenfunctions'];
        tovT=(0:sample.d)/sample.d;
        hat=@(t) 1-min(abs(t),1); %triangular hat function
        tcent=1-net(sobolset(1),sample.s); %van der Corput sequence
        xmat=zeros(sample.N,sample.d+1);       
        xmat=xmat+x(:,1)*(tovT*sqrt(sample.T));
        if sample.s>1;
            tpowm=2.^(1+floor(log2(1:sample.s-1)));
            for k=2:sample.s;
                xmat=xmat+x(:,k)*hat((tovT-tcent(k))*tpowm(k-1))*(sqrt(sample.T/(tpowm(k-1)*2)));
            end
        end
        smat=asset.s0*exp((asset.r-asset.sig.*asset.sig/2)*repmat(tovT*sample.T,sample.N,1) + xmat*asset.sig);
end


        
