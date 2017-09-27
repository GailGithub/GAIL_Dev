function plot_error_vs_time(error, time, errTol, fName, figSavePathName, ...
  errorLimits, timeLimits,visiblePlot)

% scatter Plot error vs time
if exist('visiblePlot','var') && visiblePlot==false
    hFigTime = figure('visible','off');
else
    hFigTime = figure(); 
end

%errLowLimit = 1e-7; errHighLimit = 1e1; timeLowLimit = 1e-1; timeHighLimit = 1e1;

%errLowLimit = min(error)*1e-2; errHighLimit = errTol*1e2; 
%timeLowLimit = min(time)*1e-1; timeHighLimit = max(time)*1e1;

set(hFigTime, 'units', 'inches', 'Position', [4 4 6.5 5.5])
loglog(error,time, 'b.', [errTol, errTol], [timeLimits(1) timeLimits(2)], 'r-.')

legend({'error', ...
      'errorTol'},...
      'location','northwest') 
             
axis([errorLimits(1) errorLimits(2) timeLimits(1) timeLimits(2)])
set(gca,'Xtick',(10.^(log10(errorLimits(1)):2:log10(errorLimits(2)))), ...
        'YTick',(10.^(log10(timeLimits(1)) :1:log10(timeLimits(2)))))
%set(gca,'Xtick',(10.^(-15:3:2)),'YTick',(10.^(-3:1:0)) )

xlabel('Error')
ylabel('Time (seconds)')
title(fName)

saveas(hFigTime, figSavePathName)

