function plot_error_vs_errbd(error, errBd, fName, figSavePathName, ...
  errorLimits, errbdLimits,visiblePlot)

% scatter Plot error vs time
if exist('visiblePlot','var') && visiblePlot==false
    hFigTime = figure('visible','off');
else
    hFigTime = figure(); 
end

%errLowLimit = 1e-7; errHighLimit = 1e1; timeLowLimit = 1e-1; timeHighLimit = 1e1;

set(hFigTime, 'units', 'inches', 'Position', [4 4 9.5 7.5])
loglog(error,errBd, 'b.', [errbdLimits(1), errbdLimits(2)], [errbdLimits(1) errbdLimits(2)], 'r-.')

%legend({'error', 'part'}, 'location','best') 
             
axis([errorLimits(1) errorLimits(2) errbdLimits(1) errbdLimits(2)])
set(gca,'Xtick',(10.^(log10(errorLimits(1)):3:log10(errorLimits(2)))), ...
        'YTick',(10.^(log10(errbdLimits(1)) :3:log10(errbdLimits(2)))))


xlabel('Error')
ylabel('Bayes Cub. Error Bound')
title(fName)

saveas(hFigTime, figSavePathName)

