%% Produce plots for MVNExample
gail.InitializeWorkspaceDisplay %clean up
load MVNProbExampleAllData.mat


%% IID and unscrambled Sobol
axisvec = [100 1e4 1e-12 0.1];
xtick = 10.^(2:4);
ytick = 10.^(-12:2:2);

figure

h = loglog(nvecMLE,errmedMVNProbMLESobolGn,'.', ...
      [nvecMLE(1) nvecMLE(end)],errmedMVNProbMLESobolGn(1)*[1 (nvecMLE(1)/nvecMLE(end))^2], '--', ...
      [nvecMLE nvecMLE]', [errmedMVNProbMLESobolGn errtopMVNProbMLESobolGn]' , '-', ...
      'color', noYcolorSequence{4});

legend(h([1 2]),{
   'Bayesian Cubature Sobol''', '\(O(n^{-2})\)'}, ...
   'location','best')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')

print -depsc MVNBayesianWtSobol_matern.eps


visiblePlot = false
%% Plot n vs time for bayesian cubature
if exist('visiblePlot','var') && visiblePlot==false
    hFigTime = figure('visible','off');
else
    hFigTime = figure();
end

nii = 2.^(7:13);
computeTime = [0.080798000000000, 0.104679000000000, 0.292205000000000, 1.481677000000000, ....
 14.238427000000000, 1.083312010000000e+02, 7.074992927000000e+03];

hTime = loglog(nii,computeTime, 'b.-', ...
  [nii(1) nii(end)],computeTime(1)*[1 (nii(end)/nii(1))^3 ], 'g--')
legend(hTime, {'Comp. time', '\(O(n^{3})\)'}, ...
                'location','best')
legend boxoff

% loglog(nii,computeTime, 'b.-', ...
%   nii, computeTime(1)*(exp(0.2*nii/nii(1))), 'r--', ...
%   [nii(1) nii(end)],computeTime(1)*[1 (nii(end)/nii(1))^3 ], 'g--')
% legend({'Compute time', '\(O(e^{0.2n})\)', '\(O(n^{3})\)'}, ...
%                 'location','southeast')
              
axis([100 1e4 1e-2 1e5])
set(gca,'Xtick',(10.^(2:4)),'YTick',(10.^(-2:2:5)))
xlabel('Sample Size, \(n\)')
ylabel('Comp.\ time')
%title('MVN with Matern kernel')
print -depsc MVN_bayesianCubaturecomputeTime_matern.eps
