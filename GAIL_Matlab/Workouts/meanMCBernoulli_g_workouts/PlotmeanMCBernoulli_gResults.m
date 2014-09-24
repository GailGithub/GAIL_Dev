set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
close all
clear all
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
loglog(res(maxsample,4),res(maxsample,7),'k*', ...
    'linewidth',1);
hold on;
loglog(res(standard,4),res(standard,7),'k.', ...
    'linewidth',1);
legend('Budget exceeded','Guaranteed','location','SouthEast')
hold off;

xlabel('true p')
%ylabel('abserr/abstol')
ylabel('relerr/reltol')

title('The rato of relative error and error tolerance for different p')

print('-depsc', ['./Figures/',out_param.errtype,'.eps'])
