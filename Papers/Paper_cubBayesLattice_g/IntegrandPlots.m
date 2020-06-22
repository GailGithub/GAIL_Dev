function IntegrandPlots()

% d=3 problem reduced to d=2 using Genz method
MVNParams.C = [4 1 1; 0 1 0.5; 0 0 0.25];
MVNParams.Cov = MVNParams.C'*MVNParams.C;
MVNParams.a = [-6 -2 -2];
MVNParams.b = [5 2 1];
MVNParams.mu = 0;
MVNParams.CovProp.C = chol(MVNParams.Cov)';

GenzFuncInp = @(x)GenzFunc(x,MVNParams);
%% Plot the Genz function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNParams.a) - 1;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
  
varTxVec = {'Baker', 'C0','C1','C1sin', 'C2sin', 'none' };
for varTx=varTxVec
  GenzFuncPer = cubBayesLattice_g.doPeriodTx(GenzFuncInp, varTx{1});
  zz = reshape(GenzFuncPer(xyplot),nx,nx);
  figH4 = figure();
  surf(xx,yy,zz)
  shading interp
  xlabel('\(x_1\)')
  ylabel('\(x_2\)')
  if strcmp(varTx, 'none')
    zlabel('\(f_{\textrm{Genz}}(x_1,x_2)\)')
  else
    zlabel('\(f_{\textrm{GenzP}}(x_1,x_2)\)')
  end
  view(-20,20)
  gail.save_image(figH4, 'Paper_cubBayesLattice_g', sprintf('GenzFunc_varTx_%s',varTx{1}))
end

%% Plot the Keister function
dim = 1;
  
xplot = linspace(0,1,2^10);

if dim==2
  nx = numel(xplot);
  [xx,yy] = meshgrid(xplot);
  xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
  zz = reshape(keisterFunc(xyplot,dim,0.8),nx,nx);
  figH1 = figure();
  surf(xx,yy,zz)
  shading interp
  xlabel('\(x_1\)')
  ylabel('\(x_2\)')
  zlabel('\({f}_{\textrm{Keister}}(x_1,x_2)\)')
  view(-20,20)
  gail.save_image(figH1, 'Paper_cubBayesLattice_g', sprintf('Keister_cube_a%d',ceil(a*1000)))
else
  lineShapes = {'-','--',':','-.'};

  figH1=figure(); 
  aVec = [0.8 1.5];
  for iter=1:length(aVec)
    plot(xplot', keisterFunc(xplot',1,aVec(iter)), lineShapes{iter})
    hold on
  end
  temp = string(aVec);
  temp = strcat('\(a=',temp,'\)');
  legend(temp,'location','south'); axis tight
  xlabel('\(x\)')
  ylabel('\({f}_{\textrm{Keister}}(x)\)')
  gail.save_image(figH1, 'Paper_cubBayesLattice_g', 'Keister_cube_1D')
end

% whole R domain
fKeister = @(x) cos(sum(x.^2, 2)).*exp(-sum(x.^2, 2));
if dim==2
  xplot = (-5:0.002:5);
  nx = numel(xplot);
  [xx,yy] = meshgrid(xplot);
  xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
  zz = reshape(fKeister(xyplot),nx,nx);
  figH2 = figure();
  surf(xx,yy,zz)
  shading interp
  xlabel('\(x_1\)')
  ylabel('\(x_2\)')
  zlabel('\(f_{\textrm{Keister}}(x_1,x_2)\)')
  view(-20,20)
  gail.save_image(figH2, 'Paper_cubBayesLattice_g', 'Keister_wholeR')
else
  xplot = -5:0.001:5;
  figH2 = figure();
  plot(xplot',fKeister(xplot'))
  xlabel('$x \in \bf{R}$', 'Interpreter','latex')
  ylabel('$g_{\textrm{Keister}}(x)$', 'Interpreter','latex')
  gail.save_image(figH2, 'Paper_cubBayesLattice_g', 'Keister_wholeR_1D')
end


%% Plot the Exp(Cos) function
f = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = 2;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(f(xyplot),nx,nx);
figH3 = figure();
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\( e^{\cos(2\pi x_1) + \cos(2\pi x_2)}\)')
view(-20,20)
gail.save_image(figH3, 'Paper_cubBayesLattice_g', 'ExpCos')


%% Produce plots for MVNExample

%% Plot the Genz function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNParams.a) - 1;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(GenzFuncInp(xyplot),nx,nx);
figH4 = figure();
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\(f_{\textrm{Genz}}(x_1,x_2)\)')
view(-20,20)
gail.save_image(figH4, 'Paper_cubBayesLattice_g', 'GenzFun')

end
