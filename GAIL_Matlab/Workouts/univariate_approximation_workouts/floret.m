f=@(t) sin(2 * t) .* cos(2 * t);
t=0:0.0001:2*pi;
subplot(1,2,1)
polar(t,f(t),'b')
hold on;

pause(2)

a=0;
b=2 * pi;
[q,out]=funappx_g(@(t) sin(2 * t) .* cos(2 * t), a, b);
polar(t,q(t), 'r')
norm(q(t)-f(t))
hold off

subplot(1,2,2)
err = @(x) abs(f(x) - q(x));
x=0:0.0001:2*pi;
plot(x,err(x))
axis square
axis tight