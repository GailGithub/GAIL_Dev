function [fappx,out_param,gappx,out_gparam]=testfunction(tol,nlo,nhi)
h=@(x) 5*x.^2.*(x>=0).*(x <=0.1)+(x-0.05).*(x>0.1).*(x <=0.9)...
    +(0.9-5*(x-1).^2).*(x>0.9).*(x<=1.1)+(1.95-x).*(x>1.1).*(x<=1.9)...
    +5*(x-2).^2.*(x>1.9).*(x<2);
hdoubleprime=@(x) 10*(x>=0).*(x <=0.1) + 10.*(x>0.9).*(x<=1.1)...
     +10.*(x>1.9).*(x<2);
f=@(x) h(mod(20*x,2))/400;
fdoubleprime = @(x)  hdoubleprime(mod(20*x,2));
[fappx,out_param]=funappx_g(f,0,1,tol,nlo,nhi);
t = 0:0.0001:1;
figure;
subplot(2,1,1);
plot(t,f(t));
title('test function f(x)')
subplot(2,1,2);
plot(t,abs(fdoubleprime(t)));
title('|f"(x)|')
ylim([-1,11]);
% gail.save_eps('WorkoutFunappxOutput', 'testfunction');
g=@(x) (x+8.85).*(x>=0).*(x<0.1)+(9-5*(x-0.2).^2).*(x>=0.1).*(x<0.3)...
    +(-x+9.25).*(x>=0.3).*(x<=1);
gdoubleprime=@(x) 10.*(x>=0.1).*(x<0.3);
[gappx,out_gparam]=funappx_g(g,0,1,tol,nlo,nhi);
t = 0:0.0001:1;
figure;
subplot(2,1,1);
plot(t,g(t));
title('test function g(x)')
subplot(2,1,2);
plot(t,abs(gdoubleprime(t)));
title('|g"(x)|')
ylim([-1,11]);
% gail.save_eps('WorkoutFunappxOutput', 'testfunctiong');
