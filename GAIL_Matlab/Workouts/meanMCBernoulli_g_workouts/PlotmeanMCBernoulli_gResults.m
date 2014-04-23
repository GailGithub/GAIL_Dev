set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
close all
clear all
%load TestmeanMCBernoulli-on-abs-27-Mar-2014_22.19.32.mat
%load TestmeanMCBernoulli-on-abs-27-Mar-2014_11.57.54.mat
%load TestmeanMCBernoulli-on-abs-27-Mar-2014_11.42.48.mat
%load TestmeanMCBernoulli-on-abs-23-Mar-2014_10.44.02.mat
%load TestmeanMCBernoulli-on-abs-clt-21-Mar-2014_23.30.21.mat
load TestmeanMCBernoulli-on-rel-23-Mar-2014_11.50.44.mat
%load TestmeanMCBernoulli-on-both-23-Mar-2014_15.06.33.mat
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
%legend('sample budget exceed, no guarantee')
hold on;
loglog(res(standard,4),res(standard,7),'k*', ...
    'linewidth',1);
% hold on
% loglog(res(standard,4),res(standard,1),'b');
legend('budget exceeded','guaranteed','location','SouthEast')
hold off;

xlabel('true p')
%ylabel('abserr/abstol')
ylabel('relerr/reltol')

title('The error / tolerance comparison for different p')

print('-depsc', ['./Results/' in_param.index,'.eps'])
%title('clt results')
%if res(:,5)==res(:,6)
%loglog(res(:,4),res(:,7),'r*')