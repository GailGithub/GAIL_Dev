function [x,y]=plotbatmanfunappx(abstol)
h = 0.0001;
posx = [];
posy = [];
negx = [];
negy = [];
f1= @(x) 3*sqrt(1-(x/7).^2);
a1 = [-7; 3]; b1 = [-3; 7];
[fappx,~]=funappx_g(f1,a1(1),b1(1),abstol);
t = a1(1):h:b1(1);
figure(1);
plot(t,fappx(t),'k*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
f7 = @(x) 6*sqrt(10)/7 + ...
    (1.5-0.5*abs(x)).*sqrt(abs(abs(x)-1)./(abs(x)-1)) - ...
    6*sqrt(10)/14 * sqrt(4 - (abs(x)-1).^2);
a7 = [-3 ; 1.000001]; b7 = [-1.000001; 3];
[fappx,~]=funappx_g(f7,a7(1),b7(1),abstol);
t = a7(1):h:b7(1);
plot(t,fappx(t),'m*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
f4 = @(x) 9*sqrt(abs((abs(x)-1).*(abs(x)-.75))./((1-abs(x)).*(abs(x)-.75))) - ...
    8 * abs(x);
a4 = [-0.999999;.750001]; b4=[-.750001;0.99999];
[fappx,~]=funappx_g(f4,a4(1),b4(1),abstol);
t = a4(1):h:b4(1);
plot(t,fappx(t),'g*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
f5 = @(x) 3*abs(x)+...
   0.75*sqrt(abs((abs(x)-0.75).*(abs(x)-0.5))./((0.75-abs(x)).*(abs(x)-0.5)));
a5 = [ -.749999;0.500001]; b5=[-0.500001; .749999];
[fappx,~]=funappx_g(f5,a5(1),b5(1),abstol);
t = a5(1):h:b5(1);
plot(t,fappx(t),'b*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
f6 = @(x) 2.25*sqrt(abs((x-0.5).*(x+0.5))./((0.5-x).*(0.5+x)));
a6 = -0.5; b6 = 0.5;
[fappx,~]=funappx_g(f6,a6,b6,abstol);
t = a6:h:b6;
plot(t,fappx(t),'y*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
[fappx,~]=funappx_g(f5,a5(2),b5(2),abstol);
t = a5(2):h:b5(2);
plot(t,fappx(t),'b*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
[fappx,~]=funappx_g(f4,a4(2),b4(2),abstol);
t = a4(2):h:b4(2);
plot(t,fappx(t),'g*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
[fappx,~]=funappx_g(f7,a7(2),b7(2),abstol);
t = a7(2):h:b7(2);
plot(t,fappx(t),'m*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
[fappx,~]=funappx_g(f1,a1(2),b1(2),abstol);
t = a1(2):h:b1(2);
plot(t,fappx(t),'k*');
posx = [posx t];
posy = [posy fappx(t)];
hold on;
f2 = @(x) -3*sqrt(1-(x/7).^2);
a2 = [-7; 4]; b2 = [-4; 7];
[fappx,~]=funappx_g(f2,a2(1),b2(1),abstol);
t = a2(1):h:b2(1);
plot(t,fappx(t),'k*');
negx = [negx t];
negy = [negy fappx(t)];
hold on;
f3 =@(x) abs(x/2) - ...
    ((3*sqrt(33)-7)/112)*x.^2 - ...
    3 + ...
    sqrt(1-(abs(abs(x)-2)-1).^2);
a3 = [-4;-2;0;2]; b3=[-2;0;2;4];
[fappx,~]=funappx_g(f3,a3(1),b3(1),abstol);
t = a3(1):h:b3(1);
plot(t,fappx(t),'r*');
negx = [negx t];
negy = [negy fappx(t)];
hold on;
[fappx,~]=funappx_g(f3,a3(2),b3(2),abstol);
t = a3(2):h:b3(2);
plot(t,fappx(t),'r*');
negx = [negx t];
negy = [negy fappx(t)];
hold on;
[fappx,~]=funappx_g(f3,a3(3),b3(3),abstol);
t = a3(3):h:b3(3);
plot(t,fappx(t),'r*');
negx = [negx t];
negy = [negy fappx(t)];
hold on;
[fappx,~]=funappx_g(f3,a3(4),b3(4),abstol);
t = a3(4):h:b3(4);
plot(t,fappx(t),'r*');
negx = [negx t];
negy = [negy fappx(t)];
hold on;
[fappx,~]=funappx_g(f2,a2(2),b2(2),abstol);
t = a2(2):h:b2(2);
plot(t,fappx(t),'k*');
axis equal
gail.save_eps('WorkoutFunappxOutput', 'FunAppxBatman1');

negx = [negx t];
negy = [negy fappx(t)];
figure(2)
hold on
plot(posx,posy,'LineWidth', 5);
plot(negx,negy, 'LineWidth', 5);
axis equal
gail.save_eps('WorkoutFunappxOutput', 'FunAppxBatman2');

figure(3)
rev_neg_x  = negx(end:-1:1);
rev_neg_y  = negy(end:-1:1);
x = [posx rev_neg_x];
y = [posy rev_neg_y];
fill(x,y, 'k');
axis equal
gail.save_eps('WorkoutFunappxOutput', 'FunAppxBatman3');


