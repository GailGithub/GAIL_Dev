function  [fappx,out_param]=testfunction4(tt,a,b,tol,nlo,nhi)
%tt=0.1; a=0; b=10; tol=1e-5; nlo=10; nhi=1000;
f=@(x,tt) (10*x-141*tt).*(x>=0).*(x<10*tt)...
    +(-(x-20*tt).^2/2/tt+9*tt).*(x>=10*tt).*(x<30*tt)...
    +(259*tt-10*x).*(x>=30*tt).*(x<=100*tt);
fdoubleprime=@(x,tt) -1/tt.*(x>=10*tt).*(x<30*tt);
[fappx,out_param]=funappx_g(@(x) f(x,tt),a,b,tol,nlo,nhi);
t = a:0.0001:b;
figure;
subplot(2,1,1);
plot(t,f(t,tt));
title('test function g(x)')
subplot(2,1,2);
plot(t,abs(fdoubleprime(t,tt)));
title('|g"(x)|')
ylim([-1,1.1/tt]);
gail.save_eps('WorkoutFunappxOutput', 'testfunction4');