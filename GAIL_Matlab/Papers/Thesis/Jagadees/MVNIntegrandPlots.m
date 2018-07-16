
%% Plot the Keister function
dim = 1;

if 0
  % whole R domain
  normsqd = @(t) sum(t.^2, 2); %squared l_2 norm of t
  fKeister = @(t,dim) cos( sqrt( normsqd(t) )) .* exp(-normsqd(t));
  integrand = @(x) fKeister(x,dim);
  xplot = linspace(-10,10,2^12);

elseif 0
  % domain = [0, 1]
  normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
  replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
  yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1
  fKeister = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
  integrand = @(x) fKeister(x,dim);
  
  xplot = linspace(0,1,2^12);
else
  a = 0.999;
  
  normsqd = @(t) sum(t.^2,2); %squared l_2 norm of t
  replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
  %yinv = @(t)(erfinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1

  yinv = @(t)erfinv(t);
  parta = @(t,a) a^dim *cos( a*sqrt( normsqd(t) ));
  partb = @(t,a) exp((1-a^2)*normsqd(t) );
  fKeister = @(t,dim,a) parta(t,a).*partb(t,a)*(sqrt(2*pi))^dim;
  integrand = @(x) fKeister(yinv(x),dim,a);
  
  if 0
  normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
  replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
  %yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1

  yinv = @(t)erfcinv(t);
  parta = @(nt,a) a^dim .*cos( a*sqrt( nt ));
  partb = @(nt,a) exp(nt*(1-a^2));
  fKeister = @(nt,dim,a) parta(nt,a).*partb(nt,a)*(sqrt(pi))^dim;
  integrand = @(x) fKeister(normsqd(yinv(x)),dim,a);
  end
  
  xplot = linspace(0,1,2^12);
end

if dim==2
  nx = numel(xplot);
  [xx,yy] = meshgrid(xplot);
  xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
  zz = reshape(integrand(xyplot),nx,nx);
  figure
  surf(xx,yy,zz)
  shading interp
  xlabel('\(x_1\)')
  ylabel('\(x_2\)')
  zlabel('\(f_{\textrm{Keister}}(x_1,x_2)\)')
  view(-20,20)
  print -depsc Keisteri_cube.eps
else
  figure; plot(xplot', integrand(xplot'))
end


f = @(x) cos(sum(x.^2, 2)).*exp(-sum(x.^2, 2));
xplot = (-5:0.002:5);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(f(xyplot),nx,nx);
figure
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\(f_{\textrm{Keister}}(x_1,x_2)\)')
view(-20,20)
print -depsc Keister_wholeR.eps


%% Plot the Exp(Cos) function
f = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = 2;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(f(xyplot),nx,nx);
figure
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel('\( e^{\cos(2\pi x_1) + \cos(2\pi x_2)}\)')
view(-20,20)
print -depsc ExpCos.eps


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

