format compact
g=@(x) (x+8.85).*(x>=0).*(x<0.1)+(9-5*(x-0.2).^2).*(x>=0.1).*(x<0.3)...
    +(-x+9.25).*(x>=0.3).*(x<=1);
gdoubleprime=@(x) 10.*(x>=0.1).*(x<0.3);
p1=@(x) -(x-0.5).^2 +25;
p1doubleprime=@(x) -2+0.*x;
p2=@(x) -5*(x-0.5).^2 +25;
p2doubleprime=@(x) -10+0.*x;
t = 0:0.0001:1;
tt = 0:0.0000001:1;
figure;
subplot(2,1,1);
plot(t,g(t));
title('test function g(x)')
subplot(2,1,2);
plot(t,abs(gdoubleprime(t)));
title('$|g''''(x)|$');
xlabel('$|g''''(x)|$');
ylim([-1,11]);
gail.save_eps('WorkoutFunappxOutput', 'testfunctiong');
figure;
subplot(2,1,1);
plot(t,p1(t));
title('test function \(p_1(x)\)')
subplot(2,1,2);
plot(t,abs(p1doubleprime(t)));
xlabel('$x$');
ylabel('$|p_1''''(x)|$')
title('$|p_1''''(x)|$')
ylim([0,3]);
gail.save_eps('WorkoutFunappxOutput', 'testp1');
figure;
subplot(2,1,1);
plot(t,p2(t));
title('test function \(p_2(x)\)')
subplot(2,1,2);
plot(t,abs(p2doubleprime(t)));
xlabel('$x$');
ylabel('\(|p_2''''(x)|\)')
title('\(|p_2''''(x)|\)')
ylim([0,11]);
gail.save_eps('WorkoutFunappxOutput', 'testp2');
lmax=5;
npoints = zeros(lmax,3);
error = zeros(lmax,3);
% errorsatisfy = flase(lmax,3);
nlo=10;
nhi=1000;
for l=1:lmax
   tol=10^(-5-l);
   [gappx,out_gparam]=funappxPenalty_g(g,0,1,tol,nlo,nhi);
   npoints(l,1)=out_gparam.npoints;
   error(l,1) = max(abs(g(tt)-gappx(tt)));
%    errorsatisfy(l,1) = (error(l,1)<tol);
   [pappx1,out_param1]=funappxPenalty_g(p1,0,1,tol,nlo,nhi);
   npoints(l,2)=out_param1.npoints;
   error(l,2) = max(abs(p1(tt)-pappx1(tt)));
%    errorsatisfy(l,2) = (error(l,2)<tol);
   [pappx2,out_param2]=funappxPenalty_g(p2,0,1,tol,nlo,nhi);
   npoints(l,3)=out_param2.npoints;
   error(l,3) = max(abs(p2(tt)-pappx2(tt)));
%    errorsatisfy(l,2) = (error(l,2)<tol);
end
display(npoints)

   
   