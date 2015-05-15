
a=0;
b=16*pi;

f1 = @(x) exp(.05.*x);
f2 = @(x) f1(x).*cos(x);
f3 = @(x) f1(x).*sin(x);
[q1,out1]=funappx_g(f1, a, b);
[q2,out2]=funappx_g(f2, a, b);
[q3,out3]=funappx_g(f3, a, b);

subplot(1,2,1);
x=a:0.0001:b;
plot3(f2(x),f3(x),f1(x))
axis square
axis tight

pause(2)


x=a:0.001:b;
plot3(q2(x),q3(x),q1(x),'r')
axis square
axis tight


subplot(1,2,2);
plot3(q2(x),q3(x),abs(f1(x)-q1(x))+abs(f2(x)-q2(x))+abs(f3(x)-q3(x)),'r')
axis square
axis tight