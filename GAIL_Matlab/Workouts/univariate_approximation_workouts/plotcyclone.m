function [x,y,z,err1,err2]=plotcyclone(a,b,tol)
% a = 0; b = 16*pi; tol = 1e-6; [x,y,z,err1,err2]=plotcyclone(a,b,tol);
% set error1 =|f1(x)-\hat{f1}(x)|+|f2(x)-\hat{f2}(x)|+|f3(x)-\hat{f3}(x)|
% set error2 =
% ((f1(x)-\hat{f1}(x))^2+(f2(x)-\hat{f2}(x))^2+(f3(x)-\hat{f3}(x))^2)^0.5
in_param.a = a;
in_param.b = b;
in_param.abstol = tol/3;

f1 = @(x) exp(.05.*x);
f2 = @(x) f1(x).*cos(x);
f3 = @(x) f1(x).*sin(x);
[q1,~]=funappx_g(f1, in_param);
[q2,~]=funappx_g(f2, in_param);
[q3,~]=funappx_g(f3, in_param);

figure
subplot(2,3,1);
t=a:0.001:b;
plot3(f2(t),f3(t),f1(t))
axis square
axis tight
title('Original Cyclone')

subplot(2,3,2);
% x=a:0.001:b;
x = q2(t);
y = q3(t);
z = q1(t);
plot3(x,y,z,'r')
title('Approximated Cyclone 1')
axis square
axis tight

err1 = abs(x-f2(t))+abs(y-f3(t))+abs(z-f1(t));
subplot(2,3,3);
plot(t,err1,'r')
title('Error1')
xlabel('t')
ylabel('Error1')
axis square
axis tight

in_param.abstol = sqrt(3)*tol/3;

[q1,~]=funappx_g(f1, in_param);
[q2,~]=funappx_g(f2, in_param);
[q3,~]=funappx_g(f3, in_param);

subplot(2,3,4);
t=a:0.001:b;
plot3(f2(t),f3(t),f1(t))
axis square
axis tight
title('Original Cyclone')

subplot(2,3,5);
% x=a:0.001:b;
x = q2(t);
y = q3(t);
z = q1(t);
plot3(x,y,z,'r')
title('Approximated Cyclone 2')
axis square
axis tight

err2 = sqrt((x-f2(t)).^2+(y-f3(t)).^2+(z-f1(t)).^2);
subplot(2,3,6);
plot(t,err2,'r')
title('Error2')
xlabel('t')
ylabel('Error2')
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
