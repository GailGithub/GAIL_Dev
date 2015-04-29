f = @(x) x.^2 .* sin(2*pi./x.^2);
a=.2;
b=2.5;
[q,out]=funappx_g(f, a, b);


subplot(1,2,1);
x=a:0.0001:b;
plot(x,q(x))
axis square
axis tight

subplot(1,2,2);
plot(x,abs(f(x)-q(x)))
axis square
axis tight