function  [pappx1,out_param1,pappx2,out_param2]=parabolatestfunction(tt,a,b,tol,nlo,nhi)
%tt=0.1; a=0; b=10; tol=1e-5; nlo=10; nhi=1000;
p1=@(x,tt) -(x-30*tt).^2 +(70*tt)^2;
p1doubleprime=@(x,tt) -2;
p2=@(x,tt) -1/2/tt*(x-30*tt).^2 +(70*tt)^2;
p2doubleprime=@(x,tt) -1/tt;
[pappx1,out_param1]=funappx_g(@(x) p1(x,tt),a,b,tol,nlo,nhi);
[pappx2,out_param2]=funappx_g(@(x) p2(x,tt),a,b,tol,nlo,nhi);
t = a:0.0001:b;
figure;
subplot(2,1,1);
plot(t,p1(t,tt));
title('test function p1(x)')
subplot(2,1,2);
plot(t,abs(p1doubleprime(t,tt)));
title('|p1"(x)|')
ylim([0,3]);
gail.save_eps('WorkoutFunappxOutput', 'testp1');
figure(2);
subplot(2,1,1);
plot(t,p2(t,tt));
title('test function p2(x)')
subplot(2,1,2);
plot(t,abs(p2doubleprime(t,tt)));
title('|p2"(x)|')
ylim([0,11]);
gail.save_eps('WorkoutFunappxOutput', 'testp2');
