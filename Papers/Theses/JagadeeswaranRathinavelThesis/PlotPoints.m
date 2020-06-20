%% Plot of IID and Sobol Points
gail.InitializeWorkspaceDisplay %clean up 
format long

%% Plot IID Random Points
rng(47)
n = 64;
d = 2;
tick = 0:0.25:1;
xIID = rand(n,d);
figIID = figure();
plot(xIID(:,1),xIID(:,2),'.')
xlabel('\(x_1\)')
ylabel('\(x_2\)')
% title('IID Points')
axis square
set(gca,'xtick',tick,'ytick',tick)
% print -depsc IIDPoints.eps
gail.save_image(figIID, 'JagadeeswaranRathinavelThesis', 'IIDPoints')

%% Plot Unscramled Sobol Points
xUSobol = net(sobolset(d),n);
figSobolU = figure();
plot(xUSobol(:,1),xUSobol(:,2),'.','color',MATLABOrange)
xlabel('\(x_1\)')
ylabel('\(x_2\)')
% title('Unscrambled Sobol'' Points')
axis square
set(gca,'xtick',tick,'ytick',tick)
% print -depsc USobolPoints.eps
gail.save_image(figSobolU, 'JagadeeswaranRathinavelThesis', 'USobolPoints')

%% Plot Scramled Sobol Points
xSSobol = net(scramble(sobolset(d),'MatousekAffineOwen'),n);
figSobolS = figure();
plot(xSSobol(:,1),xSSobol(:,2),'.','color',MATLABPurple)
xlabel('\(x_1\)')
ylabel('\(x_2\)')
% title('Scrambled Sobol'' Points')
axis square
set(gca,'xtick',tick,'ytick',tick)
% print -depsc SSobolPoints.eps
gail.save_image(figSobolS, 'JagadeeswaranRathinavelThesis', 'SSobolPoints')

%% Plot Lattice Points
rng(47)
xlattice = gail.lattice_gen(1,n,d);
shift = rand(1,d);
sxlat = mod(bsxfun(@plus,xlattice,shift),1); 
figH = figure();
plot(sxlat(:,1),sxlat(:,2),'.','color',MATLABGreen)
xlabel('\(x_1\)')
ylabel('\(x_2\)')
% title('Shifted Lattice Points')
axis square
set(gca,'xtick',tick,'ytick',tick)
gail.save_image(figH, 'JagadeeswaranRathinavelThesis', 'ShiftedLatticePoints')
