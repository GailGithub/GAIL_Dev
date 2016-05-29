
a=0;
b=10*pi;
f1=@(x) sin(x);
f2=@(x) cos(x);

[q1,out1]=funappxPenalty_g(f1, a, b);
[q2,out2]=funappxPenalty_g(f2, a, b);
subplot(1,2,1);
x=a:0.001:b;
plot3(q1(x),q2(x),x);
axis square


subplot(1,2,2);
plot3(q1(x),q2(x),abs(q1(x)-f1(x))+abs(q2(x)-f2(x)),'r')
axis square
axis tight