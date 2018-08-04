
%% Plot the Keister function
dim = 1;
  
if 0
  a = 0.99;
  normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
  replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
  %yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1
  yinv = @(t) erfcinv(t);

  parta = @(nt,a) a^dim .*cos( a*sqrt( nt ));
  partb = @(nt,a) exp(nt*(1-a^2));
  fKeister = @(nt,dim,a) parta(nt,a).*partb(nt,a)*(sqrt(pi))^dim;
  integrand = @(x) fKeister(normsqd(yinv(x)),dim,a);
else
  a = 0.75;
  normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
  replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
  %yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1
  yinv = @(t) gail.stdnorminv(t);

  parta = @(nt,a) a^dim .*cos( a*sqrt( nt ));
  partb = @(nt,a) exp(nt*(1-2*a^2)/2);
  fKeister = @(nt,dim,a) parta(nt,a).*partb(nt,a)*(2*pi)^(dim/2);
  integrand = @(x,a) fKeister(normsqd(yinv(x)),dim,a);

end

xplot = linspace(0,1,2^10);

if dim==2
  nx = numel(xplot);
  [xx,yy] = meshgrid(xplot);
  xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
  zz = reshape(integrand(xyplot),nx,nx);
  figH1 = figure();
  surf(xx,yy,zz)
  shading interp
  xlabel('\(x_1\)')
  ylabel('\(x_2\)')
  zlabel('\(\widehat{f}_{\textrm{Keister}}(x_1,x_2)\)')
  view(-20,20)
  saveas(figH1, sprintf('Keister_cube_a%d.png',ceil(a*1000)))
else
  pointShapes = {'o','s','d','^','v','<','>','p','h'};
  lineShapes = {'-','--',':','-.'};

  pointSize=30;
  figH1=figure(); 
  aVec = [0.75 0.9 1.2 1.9];
  for iter=1:length(aVec)
    plot(xplot', integrand(xplot',aVec(iter)), lineShapes{iter})
    hold on
  end
  temp = string(aVec);
  legend(temp,'location','best'); axis tight
  xlabel('\(x\)')
  ylabel('\(\widehat{f}_{\textrm{Keister}}(x)\)')
  saveas(figH1, sprintf('Keister_cube_1D.png'))

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
  saveas(figH2, sprintf('Keister_wholeR.png'))
else
  xplot = -5:0.001:5;
  figH2 = figure();
  plot(xplot',fKeister(xplot'))
  xlabel('$x \in \bf{R}$', 'Interpreter','latex')
  ylabel('$f_{\textrm{Keister}}(x)$', 'Interpreter','latex')
  saveas(figH2, sprintf('Keister_wholeR_1D.png'))
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
saveas(figH3, sprintf('ExpCos.png'))


%% Produce plots for MVNExample
gail.InitializeWorkspaceDisplay %clean up 
load MVNProbExampleAllData.mat
markerSize = num2cell(10*ones(size(markerSize)));

%% Plot the Genz function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNProbIIDGn.a) - 1;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(MVNProbIIDGn.f(xyplot),nx,nx);
figH4 = figure();
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\(f_{\textrm{Genz}}(x_1,x_2)\)')
view(-20,20)
saveas(figH4, sprintf('GenzFun.png'))

%% Plot the Affine function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNProbSobolAn.a);
xyplot = [xx(:) yy(:) 0.5*ones(nx*nx, dim-2)];
zz = reshape(MVNProbSobolAn.f(xyplot),nx,nx);
figH4 = figure();
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\(f_{\textrm{aff}}(x_1,x_2,1/2)\)')
view(-45,20)
saveas(figH4, sprintf('AffineFun.png'))

