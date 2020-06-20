%% Plot of IID and Sobol Points
function PlotPoints(bwcolor)
gail.InitializeDisplay %clean up
format long

if ~nargin
  bwcolor = 'color';
end

if strcmp(bwcolor,'bw')
  black = [0 0 0];
  MATLABBlue = black;
  MATLABOrange = black;
  MATLABPurple = black;
  MATLABGreen = black;
  %MATLABCyan = black;
end


%% Plot IID Random Points
rng(2947)
n = 64;
d = 2;
tick = 0:0.25:1;
xIID = rand(n,d);
figure
plot(xIID(:,1),xIID(:,2),'.','color',MATLABBlue)
xlabel('\(x_1\)')
ylabel('\(x_2\)')
title('IID Nodes')
axis square
set(gca,'xtick',tick,'ytick',tick)
%print -depsc IIDPoints.eps
gail.save_eps('MC_StoppingCriteriaOutput', 'IIDPoints');

%% Plot Unscramled Sobol Points
xUSobol = net(sobolset(d),n);
figure
plot(xUSobol(:,1),xUSobol(:,2),'.','color',MATLABOrange);
xlabel('\(x_1\)')
ylabel('\(x_2\)')
title('Unscrambled Sobol'' Nodes')
axis square
set(gca,'xtick',tick,'ytick',tick)
%print -depsc USobolPoints.eps
gail.save_eps('MC_StoppingCriteriaOutput', 'USobolPoints');

%% Plot Scramled Sobol Points
xSSobol = net(scramble(sobolset(d),'MatousekAffineOwen'),n);
figure
plot(xSSobol(:,1),xSSobol(:,2),'.','color',MATLABPurple)
xlabel('\(x_1\)')
ylabel('\(x_2\)')
title('Scrambled Sobol'' Nodes')
axis square
set(gca,'xtick',tick,'ytick',tick)
%print -depsc SSobolPoints.eps
gail.save_eps('MC_StoppingCriteriaOutput', 'SSobolPoints');

%% Plot Lattice Points
rng(47)
xlattice = gail.lattice_gen(1,n,d);
shift = rand(1,d);
sxlat = mod(bsxfun(@plus,xlattice,shift),1);
figure
plot(sxlat(:,1),sxlat(:,2),'.','color',MATLABGreen)
xlabel('\(x_1\)')
ylabel('\(x_2\)')
title('Shifted Lattice Nodes')
axis square
set(gca,'xtick',tick,'ytick',tick)
%print -depsc ShiftedLatticePoints.eps
gail.save_eps('MC_StoppingCriteriaOutput', 'ShiftedLatticePoints');

