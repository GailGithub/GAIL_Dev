%% Produce plots for MVNExample
gail.InitializeWorkspaceDisplay %clean up
load MVNProbExampleAllData.mat
markerSize = num2cell(10*ones(size(markerSize)));

%% IID and unscrambled Sobol
axisvec = [100 1e7 1e-12 0.1];
xtick = 10.^(2:7);
ytick = 10.^(-12:2:2);
figure
h = loglog(nvec,errmedMVNProbIIDGn,'.', ...
   [nvec(1) nlarge],errmedMVNProbIIDGn(1)*[1 sqrt(nvec(1)/nlarge)], '-', ...
   [nvec nvec]', [errmedMVNProbIIDGn errtopMVNProbIIDGn]', '-', ...
   'color', noYcolorSequence{1});
hold on
h = [h(1:2); loglog(nvec,errMVNProbuSobolGn,'.', ...
   [nvec(1) nlarge],errMVNProbuSobolGn(1)*[1 nvec(1)/nlarge], ':', ...
   'color', noYcolorSequence{2})];
% gail.colorMarkerLinePlot(h(1),1,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(3),2,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
legend(h,{'IID','\(O(n^{-1/2})\)', ...
   'Sobol''','\(O(n^{-1})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNIIDUSobol.eps

%% IID, unscrambled Sobol, and scrambled Sobol
figure
h = loglog(nvec,errmedMVNProbIIDGn,'.', ...
   [nvec(1) nlarge],errmedMVNProbIIDGn(1)*[1 sqrt(nvec(1)/nlarge)], '-', ...
   [nvec nvec]', [errmedMVNProbIIDGn errtopMVNProbIIDGn]', '-', ...
   'color', noYcolorSequence{1});
hold on
h = [h(1:2); loglog(nvec,errMVNProbuSobolGn,'.', ...
   [nvec(1) nlarge],errMVNProbuSobolGn(1)*[1 nvec(1)/nlarge], ':', ...
   'color', noYcolorSequence{2})];
h = [h(1:4); loglog(nvec,errmedMVNProbSobolGn,'.', ...
      [nvec(1) nlarge],errmedMVNProbSobolGn(1)*[1 (nvec(1)/nlarge)^1.5], '--', ...
      [nvec nvec]', [errmedMVNProbSobolGn errtopMVNProbSobolGn]' , '-', ...
      'color', noYcolorSequence{3})];
% gail.colorMarkerLinePlot(h(1),1,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(3),2,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(5),3,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
legend(h(1:6),{'IID','\(O(n^{-1/2})\)', ...
   'Sobol''','\(O(n^{-1})\)', ...
   'Scrambled Sobol''','\(O(n^{-3/2})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNIIDUSobolSobol.eps

%% scrambled Sobol Genz and affine
figure
h = loglog(nvec,errmedMVNProbSobolGn,'.', ...
      [nvec(1) nlarge],errmedMVNProbSobolGn(1)*[1 (nvec(1)/nlarge)^1.5], '--', ...
      [nvec nvec]', [errmedMVNProbSobolGn errtopMVNProbSobolGn]' , '-', ...
      'color', noYcolorSequence{3});
hold on
h = [h(1:2); loglog(nvec,errmedMVNProbSobolAn,'.', ...
      [nvec(1) nlarge],errmedMVNProbSobolAn(1)*[1 (nvec(1)/nlarge)^1.5], '--', ...
      [nvec nvec]', [errmedMVNProbSobolAn errtopMVNProbSobolAn]' , '-', ...
      'color', noYcolorSequence{4})];
% gail.colorMarkerLinePlot(h(1),3,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(3),4,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
legend(h([3:4 1:2]),{
   'Scr.\ Sobol'' w/ Affine', '\(O(n^{-3/2})\)' ...
   'Scr.\ Sobol'' w/ Genz','\(O(n^{-3/2})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNSobolGenzAff.eps

%% IID, unscrambled Sobol, and scrambled Sobol, and MLE cubature
figure
h = loglog(nvec,errmedMVNProbIIDGn,'.', ...
   [nvec(1) nlarge],errmedMVNProbIIDGn(1)*[1 sqrt(nvec(1)/nlarge)], '-', ...
   [nvec nvec]', [errmedMVNProbIIDGn errtopMVNProbIIDGn]', '-', ...
   'color', noYcolorSequence{1});
hold on
h = [h(1:2); loglog(nvec,errMVNProbuSobolGn,'.', ...
   [nvec(1) nlarge],errMVNProbuSobolGn(1)*[1 nvec(1)/nlarge], ':', ...
   'color', noYcolorSequence{2})];
h = [h(1:4); loglog(nvec,errmedMVNProbSobolGn,'.', ...
      [nvec(1) nlarge],errmedMVNProbSobolGn(1)*[1 (nvec(1)/nlarge)^1.5], '--', ...
      [nvec nvec]', [errmedMVNProbSobolGn errtopMVNProbSobolGn]' , '-', ...
      'color', noYcolorSequence{3})];
h = [h(1:6); loglog(nvecMLE,errmedMVNProbMLESobolGn,'.', ...
      [nvecMLE(1) nvecMLE(end)],errmedMVNProbMLESobolGn(1)*[1 (nvecMLE(1)/nvecMLE(end))^2], '--', ...
      [nvecMLE nvecMLE]', [errmedMVNProbMLESobolGn errtopMVNProbMLESobolGn]' , '-', ...
      'color', noYcolorSequence{4})];
% gail.colorMarkerLinePlot(h(1),1,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(3),2,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(5),3,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(7),4,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
legend(h([1:2:7 8]),{'IID', ...
   'Sobol''', ...
   'Scrambled Sobol''', ...
   'Bayesian Cubature Sobol''', '\(O(n^{-2})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNIIDUSobolSobolWtSobol.eps

%% Bayesian cubature only
figure

h = loglog(nvecMLE,errmedMVNProbMLESobolGn,'.', ...
      [nvecMLE(1) nvecMLE(end)],errmedMVNProbMLESobolGn(1)*[1 (nvecMLE(1)/nvecMLE(end))^2], '--', ...
      [nvecMLE nvecMLE]', [errmedMVNProbMLESobolGn errtopMVNProbMLESobolGn]' , '-', ...
      'color', noYcolorSequence{4});
% gail.colorMarkerLinePlot(h(1),1,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(3),2,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(5),3,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
% gail.colorMarkerLinePlot(h(7),4,noYcolorSequence,markerSequence,markerSize, ...
%    {'none'})
legend(h([1 2]),{
   'Bayesian Cubature Sobol''', '\(O(n^{-2})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNBayesianWtSobol.eps




%% Plot n vs time for bayesian cubature
figure
nii = 2.^(7:13);
computeTime = [0.080798000000000, 0.104679000000000, 0.292205000000000, 1.481677000000000, ....
 14.238427000000000, 1.083312010000000e+02, 7.074992927000000e+03];

loglog(nii,computeTime, 'b--')
axis([100 1e7 1e-2 1e5])
set(gca,'Xtick',(10.^(2:7)),'YTick',(10.^(-2:2:5)))
xlabel('Sample Size, \(n\)')
ylabel('Comp.\ time')
%title('MVN with Matern kernel')
print -depsc MVN_bayesianCubaturecomputeTime.eps



%% Plot Error Bounds
figure
h = loglog(errvecMVNProbMLESobolGn(:),errbdvecMBVProbMLESobolGn(:),'.', ...
   [1e-8 1e-2],[1e-8 1e-2], 'k--');
gail.colorMarkerLinePlot(h(1),4,noYcolorSequence,markerSequence,markerSize, ...
   {'none'})
axis([1e-8 1e-2 1e-8 1e-2])
xlabel('Error, \(|\mu - \hat{\mu}|\)')
ylabel('Bayesian Cub.\ Error Bound')
print -depsc MVNSobolWtSobolErrBd.eps
success = mean(errvecMVNProbMLESobolGn(:) <= errbdvecMBVProbMLESobolGn(:))

%% Plot Matern kernel
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
refpt = [0.4 0.6];
tempa = out.aMLE(end)*abs(xx(:)-refpt(1));
Kval = exp(-tempa).*(1 + tempa);
tempa = out.aMLE(end)*abs(yy(:)-refpt(2));
Kval = Kval.*exp(-tempa).*(1 + tempa);
zz = reshape(Kval,nx,nx);
figure
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel(['\(C(\textbf{x},(' num2str(refpt(1)) ',' num2str(refpt(2)) '))\)'])
view(-50,30)
print -depsc Matern.eps

%% Plot the Genz function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNProbIIDGn.a) - 1;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(MVNProbIIDGn.f(xyplot),nx,nx);
figure
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\(f_{\textrm{Genz}}(x_1,x_2)\)')
view(-20,20)
print -depsc GenzFun.eps

%% Plot the Affine function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNProbSobolAn.a);
xyplot = [xx(:) yy(:) 0.5*ones(nx*nx, dim-2)];
zz = reshape(MVNProbSobolAn.f(xyplot),nx,nx);
figure
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\(f_{\textrm{aff}}(x_1,x_2,1/2)\)')
view(-45,20)
print -depsc AffineFun.eps
