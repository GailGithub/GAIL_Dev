format compact
g=@(x) (3*x+0.25).*(x>=0).*(x<0.1)+(1-5*(x-0.4).^2).*(x>=0.1).*(x<0.7)...
    +(-3*x+2.65).*(x>=0.7).*(x<=1);
p1=@(x) -(x-0.5).^2 +25;
p2=@(x) -5*(x-0.5).^2 +25;
tt = 0:0.0000001:1;
h=@(x) 5*x.^2.*(x>=0).*(x <=0.3)+(3*x-0.45).*(x>0.3).*(x <=0.7)...
    +(2.1-5*(x-1).^2).*(x>0.7).*(x<=1.3)+(5.55-3*x).*(x>1.3).*(x<=1.7)...
    +5*(x-2).^2.*(x>1.7).*(x<2);
f=@(x) h(mod(20*x,2))/400;
lmax=5;
npoints = zeros(lmax,4);
error = zeros(lmax,4);
% errorsatisfy = flase(lmax,3);
nlo=10;
nhi=1000;
for l=1:lmax
   tol=10^(-5-l);
   [fappx,out_param]=funappxPenalty_g(f,0,1,tol,nlo,nhi);
   npoints(l,1)=out_param.npoints;
   error(l,1) = max(abs(f(tt)-fappx(tt)));
   [gappx,out_gparam]=funappxPenalty_g(g,0,1,tol,nlo,nhi);
   npoints(l,2)=out_gparam.npoints;
   error(l,2) = max(abs(g(tt)-gappx(tt)));
%    errorsatisfy(l,1) = (error(l,1)<tol);
   [pappx1,out_param1]=funappxPenalty_g(p1,0,1,tol,nlo,nhi);
   npoints(l,3)=out_param1.npoints;
   error(l,3) = max(abs(p1(tt)-pappx1(tt)));
%    errorsatisfy(l,2) = (error(l,2)<tol);
   [pappx2,out_param2]=funappxPenalty_g(p2,0,1,tol,nlo,nhi);
   npoints(l,4)=out_param2.npoints;
   error(l,4) = max(abs(p2(tt)-pappx2(tt)));
%    errorsatisfy(l,2) = (error(l,2)<tol);
end
display(npoints)