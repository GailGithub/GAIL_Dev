function PlotRatioHoeffCLT()
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
%set(0,'defaultLineMarkerSize',40) %latex axis labels
pw = -4:0.001:-1;
alpha = 10.^pw;
n =  length(alpha);
epsilon = 0.001;
N_CLTwor = ceil(norminv(1-alpha/2)./(4*epsilon^2));
N_hoeff = ceil(log(2./alpha)./(2*epsilon^2));
r_hoeffCLT = N_hoeff./N_CLTwor;
[~,~,~,MATLABVERSION] = GAILstart(false);
if usejava('jvm') || MATLABVERSION <= 7.12
    figure
    semilogx(alpha,r_hoeffCLT,'k','linewidth',2)
    xlabel('$\alpha$')
    ylabel('$n_{\mbox{Hoeff}}/n_{\mbox{CLT}}$')
    gail.save_eps('meanMCBerPaperOutput','plotHoeffCLTr');
end
end