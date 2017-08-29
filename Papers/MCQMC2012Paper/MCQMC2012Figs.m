% Make some plots for MCQMC 2012 paper
%% Garbage collection
function MCQMC2012Figs()
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
%clearvars, close all
format compact

%% Set fixed constants
alpha=0.01;
talpha=1-sqrt(1-alpha);
beta=0.01;
A1=0.3328;
A2=0.429;
%A3=25.7984;
%A3=30.54;
A3=18.1139;

%% Set eps/sigma and kappa range
tolovsigvec=-norminv(alpha/2)./sqrt(10.^(4:0.02:9)');
ntol=length(tolovsigvec);
kappavec=[2 10 100];
nkappa=length(kappavec);

%% Define functions
NCLT = @(tolovsig,alpha) ceil((norminv(alpha/2)./tolovsig).^2); %CLT sample
NCheb = @(tolovsig,alpha) ceil(1./((tolovsig.^2).*alpha)); %Chebychev sample size
%Solve for non-uniform Berry-Esseen sample size
NBEuniffunzero=@(logsqrtn,tolovsig,rho,alpha) ... %function that whose root is n
    normcdf(-exp(logsqrtn).*tolovsig) + exp(-logsqrtn)*A1*(rho+A2) ...
    - alpha/2;
NBEunif = @(tolovsig,rho,alpha) ... %function to compute n
    ceil(exp(2*fzero(@(x) NBEuniffunzero(x,tolovsig,rho,alpha),...
    log(sqrt(NCLT(tolovsig,alpha))))));
NBEnonuniffunzero=@(logsqrtn,tolovsig,rho,alpha) ...
    normcdf(-exp(logsqrtn).*tolovsig) + ...
    exp(-logsqrtn)*A3*rho./(1+(exp(logsqrtn).*tolovsig).^3) - alpha/2;
    %solve for non-uniform Berry-Esseen sample size
NBEnonunif = @(tolovsig,rho,alpha) ...
    ceil(exp(2*fzero(@(x) NBEnonuniffunzero(x,tolovsig,rho,alpha),...
    log(sqrt(NCLT(tolovsig,alpha))))));
NChebBE = @(tolovsig,rho,alpha)...
    min([NCheb(tolovsig,alpha) NBEunif(tolovsig,rho,alpha) ...
    NBEnonunif(tolovsig,rho,alpha)]);
vweight = @(alpha,beta,fudge2) ...
    sqrt(fudge2+(fudge2-1).*sqrt(alpha.*(1-beta)/(beta*(1-alpha))));
fudge2fun = @(kappa,nsigma,alpha) 1./(1 - sqrt((kappa-(nsigma-3)./(nsigma-1)).* ...
    ((1-alpha)./(alpha*nsigma))));
Nmubound = @(tolovsig,kappa,nsigma,alpha) ...
    max(nsigma,...
    NChebBE(tolovsig/vweight(alpha,beta,fudge2fun(kappa,nsigma,alpha)),...
    kappa.^(3/4),alpha));
Ntot = @(tolovsig,kappa,nsigma,alpha) ...
    nsigma+Nmubound(tolovsig,kappa,nsigma,alpha);
kurtmax = @(nsigma,alpha,fudge2) ...
    (nsigma-3)/(nsigma-1)+((alpha*nsigma)/(1-alpha))*(1-1/fudge2)^2;

%% Initialize sample sizes
NCLTvec=NCLT(tolovsigvec,alpha); %CLT sample size
NChebBEvec=zeros(ntol,nkappa); %Berry-Esseen sample size
NChebBEwhich=NChebBEvec; %Berry-Esseen sample size
Ntotvec0=NChebBEvec; %Upper bound on total sample ssize
nsigoptvec=NChebBEvec; %Optimal nsigma
Ntotoptvec=NChebBEvec; %Ntotal for optimal nsigma
fudgeoptvec=NChebBEvec; %fudge for optimal nsigma
nsigma0vec=zeros(1,nkappa);
fudge0vec=nsigma0vec;

%% Main loop
for k=1:nkappa
    
    tic
    %% Set kurtosis
    kappa=kappavec(k);
    rho=kappa^(3/4);

    %% Compute sample sizes
    nsigma0=kappa*4e3;
    nsigma0vec(k)=nsigma0;
    fudge0=sqrt(fudge2fun(kappa,nsigma0,talpha));
    fudge0vec(k)=fudge0;
    minlog10n=max(3,log10(kappa*(1-talpha)/talpha));
    for i=1:ntol
        tolovsig=tolovsigvec(i);
        [NChebBEvec(i,k), NChebBEwhich(i,k)]=NChebBE(tolovsig,rho,talpha);
        Ntotvec0(i,k)=Ntot(tolovsig,kappa,nsigma0,talpha);
        nsigopt=10.^(fminbnd(@(x) Ntot(tolovsig,kappa,10.^x,talpha),minlog10n,10));
        nsigoptvec(i,k)=round(nsigopt);
        Ntotoptvec(i,k)=Ntot(tolovsig,kappa,nsigoptvec(i,k),talpha);
    end
    fudgeoptvec(:,k)=sqrt(fudge2fun(kappa,nsigoptvec(:,k),talpha));

    toc
end

%% Plot ratios of sample sizes
[~,~,MATLABVERSION] = GAILstart(false); 
if usejava('jvm') || MATLABVERSION <= 7.12
    figure;
    h=loglog( ... %NCLTvec,NChebBEvec./repmat(NCLTvec,1,nkappa),'k-.',...
        NCLTvec,Ntotvec0./repmat(NCLTvec,1,nkappa),'k--',...
        NCLTvec,Ntotoptvec./repmat(NCLTvec,1,nkappa),'k-',...
        'linewidth',2);
    set(gca,'Linewidth',2);
    xlabel('$N_{CLT}$','Interpreter','latex','FontSize',20)
    %xlabel('${\it N}_{CLT}$')
    ylabel('Cost Ratios')
    %ylabel('{\it N_{\mu}}/{\it N}_{CLT}')
    axis([1e4 1e9 1 100])
    gail.save_eps('MCQMC2012PaperOutput/Results','MCSampleSizes');
    %print -deps 'MCSampleSizes.eps'
    
    %% Plot optima nsigma and fudge
    figure;
    line(NCLTvec,fudgeoptvec,'color','k','linestyle','--','linewidth',2)
    ax1 = gca;
    set(ax1,'XColor','k','YColor','k',...
        'XLim',[1e4 1e9],...
        'Xscale','log','Xtick',10.^(4:9),...
        'Yscale','linear','Ytick',1:0.1:1.7, ...
        'YLim',[1 1.7],'Linewidth',2);
    %xlabel('${\it N}_{CLT}$')
    xlabel('$N_{CLT}$','Interpreter','latex','FontSize',20)
    ylabel('{\it C}')
    axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k',...
        'XLim',[1e4 1e9],...
        'YLim',[1e4 1e9],...
        'LineWidth',2, ...
        'XScale','log','Yscale','log', ...
        'XTick',[],'YTick',10.^(4:9));
    ylabel('$n_{\sigma}$','Interpreter','latex','FontSize',20)
    line(NCLTvec,nsigoptvec,'color','k','linestyle','-','linewidth',2)
    %print -deps 'MCnsigmafudge.eps'
    gail.save_eps('MCQMC2012PaperOutput/Results','MCnsigmafudge');
end
end


