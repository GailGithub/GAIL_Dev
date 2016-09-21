function [x,y,z]=plotcyclone(a,b,tol)
% a = 0; b = 16*pi; tol = 1e-6; [x,y,z]=plotcyclone(a,b,tol);
% set error
% =max(max|f1(x)-\hat{f1}(x)|,max|f2(x)-\hat{f2}(x)|,max|f3(x)-\hat{f3}(x)|)

in_param.a = a;
in_param.b = b;
in_param.abstol = tol;

f1 = @(x) exp(.05.*x);
f2 = @(x) f1(x).*cos(x);
f3 = @(x) f1(x).*sin(x);
[q1,~]=funappxPenalty_g(f1, in_param);
[q2,~]=funappxPenalty_g(f2, in_param);
[q3,~]=funappxPenalty_g(f3, in_param);

figure
subplot(1,2,1);
t=a:0.001:b;
plot3(f2(t),f3(t),f1(t))
axis square
axis tight
title('Original Cyclone')

subplot(1,2,2);
% x=a:0.001:b;
x = q2(t);
y = q3(t);
z = q1(t);
plot3(x,y,z,'r')
title('Approximated Cyclone')
axis square
axis tight
gail.save_eps('WorkoutFunappxOutput', 'cyclone1');

figure(2)
subplot(3,1,1);
plot(t,abs(f2(t)-q2(t)),'r');
title('x,y,z error')
subplot(3,1,2);
plot(t,abs(f3(t)-q3(t)),'r');
subplot(3,1,3);
plot(t,abs(f1(t)-q1(t)),'r');
gail.save_eps('WorkoutFunappxOutput', 'cyclone2');
