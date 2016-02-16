function [negerr, poserr]=plotbatmanfunappx(abstol)
h = 0.001;
posx = [];
% posy = [];
negx = [];
% negy = [];
poserr = [];
negerr = [];
f1= @(x) 3*sqrt(1-(x/7).^2);
a1 = [-7; 3]; b1 = [-3; 7];
[fappx,~]=funappx_g(f1,a1(1),b1(1),abstol);
t = a1(1):h:b1(1);
figure(1);
plot(t,fappx(t),'k*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f1(t)-fappx(t))];
hold on;
% f7 = @(x) 6*sqrt(10)/7 + ...
%     (1.5-0.5*abs(x)).*sqrt(abs(abs(x)-1)./(abs(x)-1)) - ...
%     6*sqrt(10)/14 * sqrt(4 - (abs(x)-1).^2);
f7 = @(x) 6*sqrt(10)/7 + (1.5-0.5*abs(x)) -...
    6*sqrt(10)/14 * sqrt(4 - (1+x).^2);
a7 = -3; b7 = -1;
[fappx,~]=funappx_g(f7,a7,b7,abstol);
t = a7:h:b7;
plot(t,fappx(t),'m*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f7(t)-fappx(t))];
hold on;
% f4 = @(x) 9*sqrt(abs((abs(x)-1).*(abs(x)-.75))./((1-abs(x)).*(abs(x)-.75))) - ...
%     8 * abs(x);
f4 = @(x) 9 - 8 * abs(x);
a4 = [-1;.75]; b4=[-.75;1];
[fappx,~]=funappx_g(f4,a4(1),b4(1),abstol);
t = a4(1):h:b4(1);
plot(t,fappx(t),'g*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f4(t)-fappx(t))];
 hold on;
% f5 = @(x) 3*abs(x)+...
%    0.75*sqrt(abs((abs(x)-0.75).*(abs(x)-0.5))./((0.75-abs(x)).*(abs(x)-0.5)));
f5 = @(x) 3*abs(x) + 0.75;
a5 = [ -.75;0.5]; b5=[-0.5; .75];
[fappx,~]=funappx_g(f5,a5(1),b5(1),abstol);
t = a5(1):h:b5(1);
plot(t,fappx(t),'b*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f5(t)-fappx(t))];
hold on;
f6 = @(x) 2.25+0*x;
a6 = -0.5; b6 = 0.5;
[fappx,~]=funappx_g(f6,a6,b6,abstol);
t = a6:h:b6;
plot(t,fappx(t),'y*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f6(t)-fappx(t))];
hold on;
[fappx,~]=funappx_g(f5,a5(2),b5(2),abstol);
t = a5(2):h:b5(2);
plot(t,fappx(t),'b*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f5(t)-fappx(t))];
hold on;
[fappx,~]=funappx_g(f4,a4(2),b4(2),abstol);
t = a4(2):h:b4(2);
plot(t,fappx(t),'g*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f4(t)-fappx(t))];
hold on;
% f7 = @(x) 6*sqrt(10)/7 + ...
%     (1.5-0.5*abs(x)).*sqrt(abs(abs(x)-1)./(abs(x)-1)) - ...
%     6*sqrt(10)/14 * sqrt(4 - (abs(x)-1).^2);
f8 = @(x) 6*sqrt(10)/7 + (1.5-0.5*abs(x)) -...
    6*sqrt(10)/14 * sqrt(4 - (x-1).^2);
a8 = 1; b8 = 3;
[fappx,~]=funappx_g(f8,a8,b8,abstol);
t = a8:h:b8;
plot(t,fappx(t),'m*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f8(t)-fappx(t))];
hold on;
[fappx,~]=funappx_g(f1,a1(2),b1(2),abstol);
t = a1(2):h:b1(2);
plot(t,fappx(t),'k*');
posx = [posx t];
% posy = [posy fappx(t)];
poserr = [poserr abs(f1(t)-fappx(t))];
hold on;
f2 = @(x) -3*sqrt(1-(x/7).^2);
a2 = [-7; 4]; b2 = [-4; 7];
[fappx,~]=funappx_g(f2,a2(1),b2(1),abstol);
t = a2(1):h:b2(1);
plot(t,fappx(t),'k*');
negx = [negx t];
% negy = [negy fappx(t)];
negerr = [negerr abs(f2(t)-fappx(t))];
hold on;
f3 =@(x) abs(x/2) - ...
    ((3*sqrt(33)-7)/112)*x.^2 - ...
    3 + ...
    sqrt(1-(abs(abs(x)-2)-1).^2);
 a3 = -4; b3 = 4;
 [fappx,~]=funappx_g(f3,a3,b3,abstol);
  t = a3(1):h:b3(1);
 plot(t,fappx(t),'r*');
 negx = [negx t];
%  negy = [negy fappx(t)];
 negerr = [negerr abs(f3(t)-fappx(t))];
 hold on;
% a3 = [-4;-2;0;2]; b3=[-2;0;2;4];
% [fappx,~]=funappx_g(f3,a3(1),b3(1),abstol);
% t = a3(1):h:b3(1);
% plot(t,fappx(t),'r*');
% negx = [negx t];
% negy = [negy fappx(t)];
% hold on;
% [fappx,~]=funappx_g(f3,a3(2),b3(2),abstol);
% t = a3(2):h:b3(2);
% plot(t,fappx(t),'r*');
% negx = [negx t];
% negy = [negy fappx(t)];
% hold on;
% [fappx,~]=funappx_g(f3,a3(3),b3(3),abstol);
% t = a3(3):h:b3(3);
% plot(t,fappx(t),'r*');
% negx = [negx t];
% negy = [negy fappx(t)];
% hold on;
% [fappx,~]=funappx_g(f3,a3(4),b3(4),abstol);
% t = a3(4):h:b3(4);
% plot(t,fappx(t),'r*');
% negx = [negx t];
% negy = [negy fappx(t)];
% hold on;
[fappx,~]=funappx_g(f2,a2(2),b2(2),abstol);
t = a2(2):h:b2(2);
plot(t,fappx(t),'k*');
axis equal
gail.save_eps('WorkoutFunappxOutput', 'FunAppxBatman1');

negx = [negx t];
% negy = [negy fappx(t)];
negerr = [negerr abs(f2(t)-fappx(t))];
figure(2)
hold on
subplot(2,1,1)
semilogy(posx,poserr)
subplot(2,1,2);
semilogy(negx,negerr);
gail.save_eps('WorkoutFunappxOutput', 'FunAppxBatman2');

% figure(3)
% rev_neg_x  = negx(end:-1:1);
% rev_neg_y  = negy(end:-1:1);
% x = [posx rev_neg_x];
% y = [posy rev_neg_y];
% fill(x,y, 'k');
% axis equal
% gail.save_eps('WorkoutFunappxOutput', 'FunAppxBatman3');


