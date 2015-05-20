f = @(t) 2 * sin(2 * t) .* cos(2 * t);
h = 0.0001;
t = 0:h:2*pi;
subplot(1,3,1)
polar(t,f(t),'r')
hold on;
 
subplot(1,3,2)
a = 0;
b = 2 * pi;
[q,out] = funappx_g(f, a, b);
polar(t,q(t), 'b')
hold off

subplot(1,3,3)
err = @(x) abs(f(x) - q(x));
x = 0:h:2*pi;
plot(x,err(x))
axis square
axis tight
gail.save_eps('WorkoutFunappxOutput', 'WorkoutFunAppxFloret');
