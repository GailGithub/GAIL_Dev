function PlotmeanMCBer_gResults()

set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
if exist('TestmeanMCBernoulli-on-abs-06-Sep-2014_03.18.10.mat')
    load TestmeanMCBernoulli-on-abs-06-Sep-2014_03.18.10.mat
else
    warning(['TestmeanMCBernoulli-on-abs-06-Sep-2014_03.18.10.mat does not exist. '...
        'call function Test_meanMCBer_g to produce the MAT file.'])
    [~, ~, nrep, res] = Test_meanMCBer_g;
end
standard = true(nrep,1);
maxsample=res(:,5)== res(:,6);
standard = standard &~maxsample;

if usejava('jvm') || MATLABVERSION <= 7.12
    loglog(res(maxsample,4),res(maxsample,7),'k*', ...
        'linewidth',1);
    hold on;
    loglog(res(standard,4),res(standard,7),'k.', ...
        'linewidth',1);
    pos = [0.6 0.15 0.1 0.1];
    legend('BudgetExceeded','location',pos)
    hold off;
    legend boxoff
    xlabel('$p$')
    ylabel('$|p-\hat{p}_n|/\varepsilon$')
    gail.save_eps('meanMCBerPaperOutput',out_param.errtype);
end
end