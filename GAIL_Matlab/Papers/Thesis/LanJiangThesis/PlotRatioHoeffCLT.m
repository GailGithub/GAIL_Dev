close all, clear all
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
%set(0,'defaultLineMarkerSize',40) %latex axis labels
pw = -4:0.001:-1;
alpha = 10.^pw;
n =  length(alpha);
epsilon = [0.01,0.001,0.0001];
for i=1:length(epsilon)
N_CLTwor = ceil(norminv(1-alpha/2)./(4*epsilon(i)^2));
N_hoeff = ceil(log(2./alpha)./(2*epsilon(i)^2));
r_hoeffCLT(i,:) = N_hoeff./N_CLTwor;
end
figure
semilogx(alpha,r_hoeffCLT(1,:),'r','linewidth',2)
xlabel('$\alpha$')
ylabel('$N_{\mbox{Hoeff}}/N_{\mbox{CLT}}$')
print('-depsc', ['./Figures/plotHoeffCLTr','.eps'])