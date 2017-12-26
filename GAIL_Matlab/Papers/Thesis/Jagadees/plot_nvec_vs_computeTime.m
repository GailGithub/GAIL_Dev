function plot_nvec_vs_computeTime(nvec, computeTime, visiblePlot, ...
  figSavePath, samplingMethod)
% disabled for now

%end
% function dummy()

if ~exist('samplingMethod','var')
  samplingMethod = 'Lattice';
end
  
% Plot n vs time for the given cubature
if exist('visiblePlot','var') && visiblePlot==false
    hFigTime = figure('visible','off');
else
    hFigTime = figure(); 
end
nvevMaternKernel = 2.^(7:13);
computeTimeMaternKernel = [0.080798000000000, 0.104679000000000, 0.292205000000000, 1.481677000000000, ....
 14.238427000000000, 1.083312010000000e+02, 7.074992927000000e+03];

set(hFigTime, 'units', 'inches', 'Position', [1 1 6.5 5.5])
loglog(nvec,computeTime, 'b.-', ...
  nvevMaternKernel, computeTimeMaternKernel, 'r:', ...
  [nvec(1) nvec(end)],computeTime(2)*[1 (nvec(end)/nvec(1))^1], 'g--')

legend({samplingMethod, 'Matern', '\(O(n^{})\)'}, ...
                'location','best')
              
axis([100 1e6 1e-2 1e5])
set(gca,'Xtick',(10.^(2:6)),'YTick',(10.^(-2:2:5)))
xlabel('Sample Size, \(n\)')
ylabel('Comp.\ time')
%title('MVN with Matern kernel')

saveas(hFigTime, figSavePath)

%plot_nvec_vs_computeTime_matern_example(visiblePlot)


end

