set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
close all
clear all
%load TestmeanMCBernoulli-on-abs-03-Dec-2015_18.50.43.mat
load TestmeanMCBernoulli-on-rel-09-Sep-2014_00.15.30.mat
%load TestmeanMCBernoulli-on-abs-06-Sep-2014_03.18.10.mat
standard = true(nrep,1);
maxsample=res(:,5)== res(:,6);
%maxsample=res(:,7)== res(:,8);
standard = standard &~maxsample;
% loglog(res(maxsample,5),res(maxsample,9),'k*', ...
%     res(maxsample,5),res(maxsample,10),'go',...
%     'linewidth',1);
% hold on;
% loglog(res(standard,5),res(standard,9),'r*', ...
%     res(standard,5),res(standard,10),'bo',...
%     'linewidth',1);
% hold off;
loglog(res(maxsample,4),res(maxsample,7),'r*', ...
    'linewidth',1);
hold on;
loglog(res(standard,4),res(standard,7),'k.', ...
    'linewidth',1);
pos = [0.6 0.15 0.1 0.1];
legend('BudgetExceeded','location',pos)
hold off;
legend boxoff
xlabel('$p$')
ylabel('$|p-\hat{p}_n|/(p\varepsilon_r)$')



print('-depsc', ['./Figures/',out_param.errtype,'.eps'])
