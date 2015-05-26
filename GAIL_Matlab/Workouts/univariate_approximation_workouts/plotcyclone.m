function [x,y,z,err]=plotcyclone(a,b,tol)
% a = 0; b = 16*pi; tol = 1e-6; [x,y,z,err]=plotcyclone(a,b,tol);
in_param.a = a;
in_param.b = b;
in_param.abstol = tol;

f1 = @(x) exp(.05.*x);
f2 = @(x) f1(x).*cos(x);
f3 = @(x) f1(x).*sin(x);
[q1,~]=funappx_g(f1, in_param);
[q2,~]=funappx_g(f2, in_param);
[q3,~]=funappx_g(f3, in_param);

figure
subplot(1,3,1);
t=a:0.001:b;
plot3(f2(t),f3(t),f1(t))
axis square
axis tight
title('Original Cyclone')

subplot(1,3,2);
% x=a:0.001:b;
x = q2(t);
y = q3(t);
z = q1(t);
plot3(x,y,z,'r')
title('Approximated Cyclone')
axis square
axis tight

err = sqrt((x-f2(t)).^2+(y-f3(t)).^2+(z-f1(t)).^2);
subplot(1,3,3);
plot(t,err,'r')
title('Error')
xlabel('t')
ylabel('MSE')
axis square
axis tight
gail.save_eps('WorkoutFunappxOutput', 'WorkoutFunAppxCyclone');

% figure(2)
% subplot(1,3,1);
% plot(x,abs(f2(x)-q2(x)),'r');
% subplot(1,3,2);
% plot(x,abs(f3(x)-q3(x)),'r');
% subplot(1,3,3);
% plot(x,abs(f1(x)-q1(x)),'r');
