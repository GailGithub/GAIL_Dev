function WalshFourierCoeffDecay
[~,~,MATLABVERSION] = GAILstart(false);

if usejava('jvm') || MATLABVERSION <= 7.12
    %% Garbage collection and initialization
    format compact %remove blank lines from output
    format long e %lots of digits
    clear all %clear all variables
    close all %close all figures
    set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
    set(0,'defaultLineLineWidth',3) %thick lines
    set(0,'defaultTextInterpreter','latex') %latex axis labels
    set(0,'defaultLineMarkerSize',40) %latex axis labels
    printc = 'color'; % choose between color and bnw

    %% Initialize parameters
    mmax=20; %maximum number of points is 2^mmax
    mdualvec=12:12;
    mplot=16;
    % mmax=16; %maximum number of points is 2^mmax
    % mdualvec=11;
    % mplot=14;
    mlag=4;
    testfun=@(x) exp(-3*x).*sin(10*x.^2); d=1; %test function
    sobstr=sobolset(d);
    sobstr=scramble(sobstr,'MatousekAffineOwen');
    sobol=qrandstream(sobstr);

    %% Plot function
    figure
    xplot=(0:0.002:1);
    yplot=testfun(xplot);
    plot(xplot,yplot,'b-');
    ymin=1.1*min(yplot);
    ymax=1.1*max(yplot);
    axis([0 1 ymin ymax])
    gail.save_eps('Paper_cubSobol_g', 'Paper_cubSobol_g_FunctionWalshFourierCoeffDecay');

    %% Evaluate Function and FWT
    n=2^mmax;
    xpts=sobstr(1:n,1:d);
    y=testfun(xpts);
    yval=y;
    %yfwt=fwht(y);

    %% Compute initial FWT
    for l=0:mmax-1
       nl=2^l;
       nmmaxlm1=2^(mmax-l-1);
       ptind=repmat([true(nl,1); false(nl,1)],nmmaxlm1,1);
       evenval=y(ptind);
       oddval=y(~ptind);
       y(ptind)=(evenval+oddval)/2;
       y(~ptind)=(evenval-oddval)/2;
    end

    %% Create kappanumap
    kappanumap=(1:n)'; %initialize map
    for l=mmax-1:-1:1
       nl=2^l;
       oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
       newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
       flip=find(newone>oldone); %
       temp=kappanumap(nl+1+flip);
       kappanumap(nl+1+flip)=kappanumap(1+flip);
       kappanumap(1+flip)=temp;
       %keyboard
    end
    ymap=y(kappanumap);


    %% Plot FW coefficients
    ltgray=0.8*ones(1,3);
    gray=0.5*ones(1,3);
    nplot=2^mplot;
    yfwtabs=abs(ymap(1:nplot));
    ymin=max(1e-15,min(yfwtabs));
    ymax=max([1; yfwtabs]);
    for mdual=mdualvec
       ndual=2^mdual;
       whdual=ndual*(1:2^(mplot-mdual)-1);
       whsmall=1:ndual-1;
       whbig=ndual:nplot-1;
       muse=mdual-mlag;
       nuse=2^muse;
       whuse=nuse/2:nuse-1;
       figure
       switch printc
           case 'color'
               h=loglog(whsmall,yfwtabs(whsmall+1),'g.',...
                  whbig,yfwtabs(whbig+1),'k.',...
                  whuse,yfwtabs(whuse+1),'b.',...
                  whdual,yfwtabs(whdual+1),'r.','MarkerSize',10);
               set(h([3 4]),'MarkerSize',20)
           case 'bnw'
               h=zeros(4,1);
               h(1)=loglog(whsmall,yfwtabs(whsmall+1),'.',...
                  'MarkerSize',10,'MarkerFaceColor',ltgray,'MarkerEdgeColor',ltgray);
               hold on
               h(2)=loglog(whbig,yfwtabs(whbig+1),'.','MarkerSize',10,...
                  'MarkerFaceColor',gray,'MarkerEdgeColor',gray);
               h(3)=loglog(whuse,yfwtabs(whuse+1),'sk','MarkerSize',7,...
                  'MarkerFaceColor','k');
               h(4)=loglog(whdual,yfwtabs(whdual+1),'.k','MarkerSize',30);
       end
       maxexp=floor(log10(nplot-1));
       set(gca,'Xtick',10.^(0:maxexp))
       axis([1 nplot-1 ymin ymax])
       xlabel('$\kappa$')
       ylabel('$|\hat{f}_{\kappa}|$','interpreter','latex')
       legend(h([4 2 3]),{['error $\le \hat{S}_{0,' int2str(mdual) '}(f)$'],...
          ['$\check{S}_{' int2str(mdual) '}(f)$'],...
          ['$S_{' int2str(mdual-mlag) '}(f)$']},...
          'location','southwest','interpreter','latex')
       legend('boxoff')
       set(gca,'Position',[0.2 0.155 0.75 0.77])
       gail.save_eps('Paper_cubSobol_g', 'Paper_cubSobol_g_WalshFourierCoeffDecay');
    end
    close all
end
end
