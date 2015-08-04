function [fappx,out_param,gappx,out_gparam]=testfunctionnew(tol,nlo,nhi)
h=@(x) 5*x.^2.*(x>=0).*(x <=0.2)+(2*x-0.2).*(x>0.2).*(x <=0.8)...
    +(1.6-5*(x-1).^2).*(x>0.8).*(x<=1.2)+(3.8-2*x).*(x>1.2).*(x<=1.8)...
    +5*(x-2).^2.*(x>1.8).*(x<2);
f=@(x) h(mod(20*x,2))/400;
figure;
xplot=0:0.0002:1;
fplot=f(xplot);
plot(xplot,fplot)
figure
xpplot=xplot(1:end-1);
fpplot=diff(fplot)./diff(xplot);
plot(xpplot,fpplot)
figure
fppplot=diff(fpplot)./diff(xpplot);
plot(xpplot(1:end-1),fppplot)
%[fappx,out_param]=funappx_g(f,0,1,tol,nlo,nhi);
[fappx,out_param]=funappx_g(f,0,1,tol,nlo,nhi);
g=@(x) (2*x+0.4).*(x>=0).*(x<0.2)+(1-5*(x-0.4).^2).*(x>=0.2).*(x<0.6)...
    +(-2*x+2).*(x>=0.6).*(x<=1);
figure;
xplot=0:0.0002:1;
gplot=g(xplot);
plot(xplot,gplot)
figure
xpplot=xplot(1:end-1);
gpplot=diff(gplot)./diff(xplot);
plot(xpplot,gpplot)
figure
gppplot=diff(gpplot)./diff(xpplot);
plot(xpplot(1:end-1),gppplot)
[gappx,out_gparam]=funappx_g(g,0,1,tol,nlo,nhi);
