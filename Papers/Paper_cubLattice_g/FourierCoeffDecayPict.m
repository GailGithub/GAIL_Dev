function FourierCoeffDecayPict

%% Garbage collection and initialization

format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',15,'defaulttextfontsize',16) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels
printc = 'color'; % choose between color and bnw

%% Initialize parameters
mmax=17; %maximum number of points is 2^mmax
mdualvec=12:12; %plots of coefficients S(0,mdualvec)
mplot=16;
% mmax=16; %maximum number of points is 2^mmax
% mdualvec=11;
% mplot=14;
mlag=4;
testfun=@(x) sin(20*sqrt((x(:,1)-0.5).^2+(x(:,2)-0.5).^2)+eps)./(20*sqrt((x(:,1)-0.5).^2+(x(:,2)-0.5).^2)+eps); d=2; %test function


%% Plot the function
[~,~,MATLABVERSION] = GAILstart(false);
if usejava('jvm') || MATLABVERSION <= 7.12
[X,Y] = meshgrid(0:.01:1);
R = sqrt((X-0.5).^2 + (Y-0.5).^2) + eps;
Z = sin(20*R)./(20*R);
zmin=1.1*min(min(Z));
zmax=1.1*max(max(Z));
axis([0 1 0 1 zmin zmax])
mesh(X,Y,Z)
title('$\frac{\sin\left(20\sqrt{x^2+y^2}\right)}{20\sqrt{x^2+y^2}}$')
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_FunctionFourierCoeffDecay');
end

shift=rand;

%% Evaluate Function and FFT
n=2^mmax;
xpts=mod(gail.lattice_gen(1,n,d)+shift,1);
y=testfun(xpts);

%% Compute initial FFT
for l=0:mmax-1
   nl=2^l;
   nmmaxlm1=2^(mmax-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmmaxlm1,1);
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
   coefv=repmat(coef,nmmaxlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;
end

%% Create kappanumap
kappanumap=(1:n)'; %initialize map
for l=mmax-1:-1:1
    nl=2^l;
    oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
    newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
    flip=find(newone>oldone);
    temp=kappanumap(nl+1+flip);
    kappanumap(nl+1+flip)=kappanumap(1+flip);
    kappanumap(1+flip)=temp;
end
ymap=y(kappanumap);


%% Plot FW coefficients
[~,~,MATLABVERSION] = GAILstart(false);
if usejava('jvm') || MATLABVERSION <= 7.12
    ltgray=0.8*ones(1,3);
    gray=0.5*ones(1,3);
    nplot=2^mplot;
    yfftabs=abs(ymap(1:nplot));
    ymin=max(1e-15,min(yfftabs));
    ymax=max([1; yfftabs]);
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
                h=loglog(whsmall,yfftabs(whsmall+1),'g.',...
                    whbig,yfftabs(whbig+1),'k.',...
                    whuse,yfftabs(whuse+1),'b.',...
                    whdual,yfftabs(whdual+1),'r.','MarkerSize',10);
                set(h([3 4]),'MarkerSize',20)
            case 'bnw'
                h=zeros(4,1);
                h(1)=loglog(whsmall,yfftabs(whsmall+1),'.',...
                    'MarkerSize',10,'MarkerFaceColor',ltgray,'MarkerEdgeColor',ltgray);
                hold on
                h(2)=loglog(whbig,yfftabs(whbig+1),'.','MarkerSize',10,...
                    'MarkerFaceColor',gray,'MarkerEdgeColor',gray);
                h(3)=loglog(whuse,yfftabs(whuse+1),'sk','MarkerSize',7,...
                    'MarkerFaceColor','k');
                h(4)=loglog(whdual,yfftabs(whdual+1),'.k','MarkerSize',30);
        end
        maxexp=floor(log10(nplot-1));
        set(gca,'Xtick',10.^(0:maxexp))
        axis([1 nplot-1 ymin ymax])
        xlabel('$\kappa$')
        ylabel('$|\hat{f}_{\kappa}|$','interpreter','latex')
        legend(h([4 2 3]),{['error $\le \hat{S}_{0,' int2str(mdual) '}(f)$'],...
            ['$\check{S}_{' int2str(mdual) '}(f)$'],...
            ['$S_{' int2str(mdual-mlag) '}(f)$']},...
            'location','southwest' ,'interpreter','latex')
        legend('boxoff')
        set(gca,'Position',[0.2 0.155 0.75 0.67])
        title({'Fourier coefficients decay for $\frac{\sin\left(20\sqrt{x^2+y^2}\right)}{20\sqrt{x^2+y^2}}$'})
        gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_FourierCoeffDecay');
    end
end
close all
