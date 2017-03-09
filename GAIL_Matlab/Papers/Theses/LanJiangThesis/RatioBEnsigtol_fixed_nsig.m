% RatioBEnsigtol_fixed_nsig produces Figure 3.1 in Lan Jiang's thesis

clear all; close all;clc;
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',20) %latex axis labels
alpha = 0.01;
%tol = 0.01;
pw = 1:0.1:5;
sigtol = 10.^pw;
sigtol = sigtol';
sigtol2=sigtol.^2;
alpha1 = 1-sqrt(1-alpha);
nsig = 1e3;
A=18.1139;
A1=0.3322;
A2=0.429;
A3= 0.3031;
A4= 0.646;
A5= 0.469;
fudge = [1.01 1.1 1.2 1.3 2];
N_CLT = ceil((gail.stdnorminv(1-alpha1/2)*sigtol).^2);
ncheb = ceil(1/alpha1*sigtol.^2);
kurtmax = zeros(1,length(fudge));
logsqrtnCLT = zeros(length(sigtol),length(fudge));
logsqrtnBE = zeros(length(sigtol),length(fudge));
for j=1:length(fudge)
kurtmax(j) = kurtosismax(nsig,alpha1,fudge(j));
M3upper=kurtmax(j)^(3/4);
for i = 1:length(sigtol)
logsqrtnCLT(i,j)=log(gail.stdnorminv(1-alpha1/2).*sigtol(i));
BEfun=@(logsqrtn)gail.stdnormcdf(-exp(logsqrtn)./sigtol(i))...
    +exp(-logsqrtn).*min([A1*(M3upper+A2),A3*(M3upper+A4)...
    ,A5*M3upper, ...
    A*M3upper./(1+(exp(logsqrtn)./sigtol(i)).^3)])- alpha1/2;
logsqrtnBE(i,j) = fzero(BEfun,logsqrtnCLT(i,j));
N_BE(i,j)=ceil(exp(2*logsqrtnBE(i,j)));
%N_BE(i,j)=min(ncheb(i),ceil(exp(2*logsqrtnBE(i,j))));
end
end
[~,~,MATLABVERSION] = GAILstart(false);if usejava('jvm') || MATLABVERSION <= 7.12
    figure
        loglog(sigtol2,(N_BE(:,1)+nsig)./sigtol2,'r',sigtol2,...
            (N_BE(:,2)+nsig)./sigtol2,'k',...
        sigtol2,(N_BE(:,3)+nsig)./sigtol2,'b',sigtol2,...
        (N_BE(:,4)+nsig)./sigtol2,'g',...
        sigtol2,(N_BE(:,5)+nsig)./sigtol2,'m','linewidth',2)
    %loglog(N_CLT,fudge11,'b',N_CLT,fudge12,'g',N_CLT,fudge13,'r','linewidth',2)
    xlabel('$(\sigma/\varepsilon_a)^2$')
    ylabel('$N_{\mbox{BETotal}}$ /$(\sigma/\varepsilon_a)^2$')
    axis([1e3 1e10 1 1e4 ])
    ax = gca;
     set(ax,'XTick',[1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 ])
      set(ax,'YTick',[1 10 100 1e3 1e4])
 %set(ax,'YTick',[1 10 20 30 40 50 60 70 80 90 100])
% set(ax,'XTick',[1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 ])
% set(ax,'YTick',[ 1 2 3 4 5 6 7 8])
legend({['SDIF = ' num2str(fudge(1))],...
    ['SDIF = ' num2str(fudge(2))],...
    ['SDIF = ' num2str(fudge(3))],['SDIF = ' num2str(fudge(4))],...
    ['SDIF = ' num2str(fudge(5))]},'Location','northeast')
%legend('SDIF = 1.1','SDIF = 1.2','SDIF = 1.3')
    gail.save_eps('LanThesisOutput',['Ratio_BEn_sigtol_fixed_nsig',num2str(nsig)]);
end