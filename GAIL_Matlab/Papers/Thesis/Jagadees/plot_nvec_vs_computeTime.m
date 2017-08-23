function plot_nvec_vs_computeTime(nvec, computeTime, visiblePlot, figSavePath)
% disabled for now
end

function dummy()
  
% Plot n vs time for the given cubature
if exist('visiblePlot','var') && visiblePlot==false
    hFigTime = figure('visible','off');
else
    hFigTime = figure(); 
end

set(hFigTime, 'units', 'inches', 'Position', [4 4 6.5 5.5])
loglog(nvec,computeTime, 'b.-', ...
  [nvec(1) nvec(end)],computeTime(2)*[1 (nvec(end)/nvec(1))^1], 'g--')

legend({'Compute time', '\(O(n^{})\)'}, ...
                'location','best')
              
axis([100 1e7 1e-2 1e5])
set(gca,'Xtick',(10.^(2:7)),'YTick',(10.^(-2:2:5)))
xlabel('Sample Size, \(n\)')
ylabel('Comp.\ time')
%title('MVN with Matern kernel')

saveas(hFigTime, figSavePath)




end


function plot_nvec_vs_computeTime_matern_example(visiblePlot)

%% Plot n vs time for bayesian cubature
if exist('visiblePlot','var') && visiblePlot==false
    hFigTime = figure('visible','off');
else
    hFigTime = figure(); %21
end

nii = 2.^(7:13);
computeTime = [0.080798000000000, 0.104679000000000, 0.292205000000000, 1.481677000000000, ....
 14.238427000000000, 1.083312010000000e+02, 7.074992927000000e+03];

loglog(nii,computeTime, 'b.-', ...
  nii, computeTime(1)*(exp(0.2*nii/nii(1))), 'r--', ...
  [nii(1) nii(end)],computeTime(1)*[1 (nii(end)/nii(1))^3 ], 'g--')
legend({'Compute time', '\(O(e^{0.2n})\)', '\(O(n^{3})\)'}, ...
                'location','southeast')
              
axis([100 1e7 1e-2 1e5])
set(gca,'Xtick',(10.^(2:7)),'YTick',(10.^(-2:2:5)))
xlabel('Sample Size, \(n\)')
ylabel('Comp.\ time')
%title('MVN with Matern kernel')
print -depsc MVN_bayesianCubaturecomputeTime.eps
end
