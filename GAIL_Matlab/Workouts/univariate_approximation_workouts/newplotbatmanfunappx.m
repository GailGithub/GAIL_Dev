function [negerr, poserr]=newplotbatmanfunappx(abstol)
f= @(x) 3*sqrt(1-(x/7).^2).*(abs(x)>3).*(abs(x)<=7)+ (6*sqrt(10)/7 +...
    (1.5-0.5*abs(x)) - 6*sqrt(10)/14 * sqrt(4 - (abs(x)-1).^2))...
    .*(abs(x)<=3).*(abs(x)>1)+(9 - 8 * abs(x)).*(abs(x)<=1).*(abs(x)>0.75)...
    + (3*abs(x) + 0.75).*(abs(x)<=0.75).*(abs(x)>0.5)+(2.25+0*x).*(abs(x)<=0.5);
f2 = @(x) -3*sqrt(1-(x/7).^2).*(abs(x)<=7).*(abs(x)>4)+...
    (abs(x/2) - ((3*sqrt(33)-7)/112)*x.^2 - 3 + sqrt(1-(abs(abs(x)-2)-1).^2))...
    .*(abs(x)<=4);
[fappx,~]=funappxPenalty_g(f,-7,7,abstol,10,15);
[f2appx,~]=funappxPenalty_g(f2,-7,7,abstol,10,15);
t = -7:0.001:7;
figure(1)
plot(t,fappx(t),'b*',t,f2appx(t),'b*');
axis equal
gail.save_eps('WorkoutFunappxOutput', 'NewFunAppxBatman1');
poserr = abs(fappx(t)-f(t));
negerr = abs(f2appx(t)-f2(t));
figure(2)
subplot(2,1,1)
semilogy(t,poserr,'b-',t,ones(size(t))*abstol,'r--')
subplot(2,1,2);
semilogy(t,negerr,'b-',t,ones(size(t))*abstol,'r--');
gail.save_eps('WorkoutFunappxOutput', 'NewFunAppxBatman2');
