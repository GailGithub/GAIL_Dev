% RatioBEnsigtol_fixed_kmax produce Figure 3.2 in Lan Jiang's thesis

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
alpha1 = 1-sqrt(1-alpha);
%nsig = 1e4;
A=18.1139;
A1=0.3322;
A2=0.429;
A3= 0.3031;
A4= 0.646;
A5= 0.469;
fudge = [1.01 1.1 1.2 1.3 2];
kmax = 2;
M3upper=kmax^(3/4);
N_CLT = ceil((gail.stdnorminv(1-alpha1/2)*sigtol).^2);
logsqrtnCLT = zeros(length(sigtol),length(fudge));
logsqrtnBE = zeros(length(sigtol),length(fudge));
for j=1:length(fudge)
nsig(j) = kurtosismaxinv(kmax,alpha1,fudge(j));
%kurtmax(j) = kurtosismax(nsig,alpha1,fudge(j));
for i = 1:length(sigtol)
logsqrtnCLT(i,j)=log(gail.stdnorminv(1-alpha1/2).*sigtol(i));
BEfun=@(logsqrtn)gail.stdnormcdf(-exp(logsqrtn)./(sigtol(i)*fudge(j)))...
    +exp(-logsqrtn).*min([A1*(M3upper+A2),A3*(M3upper+A4)...
    ,A5*M3upper, ...
    A*M3upper./(1+(exp(logsqrtn)./(sigtol(i)*fudge(j))).^3)])- alpha1/2;
logsqrtnBE(i,j) = fzero(BEfun,logsqrtnCLT(i,j));
N_BE(i,j)=ceil(exp(2*logsqrtnBE(i,j)));
N_BE_tot (i,j) = N_BE(i,j)+nsig(j);
end
end
[~,~,MATLABVERSION] = GAILstart(false);
if usejava('jvm') || MATLABVERSION <= 7.12
    figure
    sigtol2=sigtol.^2;
    loglog(sigtol2,N_BE_tot(:,1)./sigtol2,'r',sigtol2,N_BE_tot(:,2)./sigtol2,'k',...
        sigtol2,N_BE_tot(:,3)./sigtol2,'b',sigtol2,N_BE_tot(:,4)./sigtol2,'g',...
        sigtol2,N_BE_tot(:,5)./sigtol2,'m','linewidth',2)
    xlabel('$(\sigma/\varepsilon_a)^2$')
    ylabel('$N_{\mbox{BETotal}}$ /$(\sigma/\varepsilon_a)^2$')
    axis([1e3 1e10 1 1e4 ])
    ax = gca;
 set(ax,'XTick',[1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 ])
 set(ax,'YTick',[1 10 100 1e3 1e4])
legend({['SDIF = ' num2str(fudge(1))],...
    ['SDIF = ' num2str(fudge(2))],...
    ['SDIF = ' num2str(fudge(3))],['SDIF = ' num2str(fudge(4))],...
    ['SDIF = ' num2str(fudge(5))]},'Location','northeast')
    gail.save_eps('LanThesisOutput',['Ratio_BEn_sigtol_fixed_kmax_',num2str(kmax)]);
end



