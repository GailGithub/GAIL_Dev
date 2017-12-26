
%% Plot the Keister function
dim = 2;
normsqd = @(t) sum(t.^2, 2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1
fKeister = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi)/2)^dim;
f = @(x) fKeister(x,dim);

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
print -depsc Keisteri_cube.eps


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

