
close all

f=@(x) floor(x/2).*(x - floor(x/2) - 1/2) + ...
   (1/2)*(mod(x,2)-1).^2 .* (mod(x,2)>=1);
xplot=0:0.002:10;
fplot=f(xplot);
plot(xplot,fplot)
figure
xpplot=xplot(1:end-1);
fpplot=diff(fplot)./diff(xplot);
plot(xpplot,fpplot)
figure
fppplot=diff(fpplot)./diff(xpplot);
plot(xpplot(1:end-1),fppplot)

a=40;
b=100;
nlo=10;
nhi=10;
tol=1e-2;
f2 =  @(x) b*f(a*x)/(a.*a);
funappx_g_gui(f2,0,1,tol,nlo,nhi)
xplot=0:0.002:1;
figure
f2plot=f2(xplot);
plot(xplot,f2plot)
figure
xpplot=xplot(1:end-1);
f2pplot=diff(f2plot)./diff(xplot);
plot(xpplot,f2pplot)
figure
f2ppplot=diff(f2pplot)./diff(xpplot);
plot(xpplot(1:end-1),f2ppplot)
[f2appx,out_paramf]=funappx_g(f2,0,1,tol,nlo,nhi);
out_paramf.npoints

g2 = @(x)  b*f(2*x)/4;
funappx_g_gui(g2,0,1,tol,nlo,nhi)
[g2appx,out_paramg]=funappx_g(g2,0,1,tol,nlo,nhi);
out_paramg.npoints
figure
g2plot=g2(xplot);
plot(xplot,g2plot)
figure
g2pplot=diff(g2plot)./diff(xplot);
plot(xpplot,g2pplot)
figure
g2ppplot=diff(g2pplot)./diff(xpplot);
plot(xpplot(1:end-1),g2ppplot)

