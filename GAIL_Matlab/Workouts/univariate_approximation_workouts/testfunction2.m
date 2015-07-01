function  [fappx,out_param]=testfunction2(tt,a,b,tol,nlo,nhi)
%tt=0.1; a=0; b=10; tol=1e-5; nlo=10; nhi=1000;
f=@(x,epsilon) (mod(x,20*epsilon).^2/2/epsilon).*(mod(x,20*epsilon)>=0).*(mod(x,20*epsilon)<epsilon)...
    +(mod(x,20*epsilon)-epsilon/2).*(mod(x,20*epsilon)>=epsilon).*(mod(x,20*epsilon)<9*epsilon)...
    +(-(mod(x,20*epsilon)-10*epsilon).^2/2/epsilon+9*epsilon).*(mod(x,20*epsilon)>=9*epsilon).*(mod(x,20*epsilon)<11*epsilon)...
    +(39*epsilon/2-mod(x,20*epsilon)).*(mod(x,20*epsilon)>=11*epsilon).*(mod(x,20*epsilon)<19*epsilon)...
    +((mod(x,20*epsilon)-20*epsilon).^2/2/epsilon).*(mod(x,20*epsilon)>=19*epsilon).*(mod(x,20*epsilon)<20*epsilon);
fdoubleprime=@(x,epsilon) 1/epsilon.*(mod(x,20*epsilon)>=0).*(mod(x,20*epsilon)<epsilon)...
    +(-1/epsilon).*(mod(x,20*epsilon)>=9*epsilon).*(mod(x,20*epsilon)<11*epsilon)...
    +1/epsilon.*(mod(x,20*epsilon)>=19*epsilon).*(mod(x,20*epsilon)<20*epsilon);
t = a:0.0001:b;
[fappx,out_param]=funappx_g(@(x) f(x,tt),a,b,tol,nlo,nhi);
figure;
subplot(2,1,1);
plot(t,f(t,tt));
title('test function f(x)')
subplot(2,1,2);
plot(t,abs(fdoubleprime(t,tt)));
title('|f"(x)|')
ylim([-1,1.1/tt]);
gail.save_eps('WorkoutFunappxOutput', 'testfunction2');
