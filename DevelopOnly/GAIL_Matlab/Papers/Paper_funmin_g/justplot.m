% Just a sample not based on the real data

close all
%% initial curve
x=[-3,-2,-1,0,1,2,3];
y=[1.2,0.4,1,0,0.3,0.5,0.1];
xx = -3:1/1000:3;
yy = spline(x,y,xx);
yy1=interp1(x,y,xx);

%% create lower bound functions
gn1=@(x) 2.4*x.^2+11.2.*x+13.2; 
gn2=@(x) 2.4*x.^2+7.8.*x+6.4; 
gn3=@(x) 5/3*x.^2+2/3.*x; 
gn4=@(x) 0.9*x.^2-0.6.*x; 
gn5=@(x) 1.2*x.^2-3.4.*x+2.5; 
gn6=@(x) 4/3*x.^2-106/15*x+9.3; 
minpt_x=[-7/3 -13/8 -1/5   1/3  17/12  53/20];
minpt_y=[2/15 1/16  -1/15 -0.1  1.1/12 -19/300];

%% plot
plot(xx,0,'black','LineWidth',3); 
h_text=text(xx(end),0,'$U_n$');
set(h_text,'FontSize',30,'interpreter','latex');
hold on
% An=plot(xx,yy1,'red');
% An_legend=legend(An,'An(f)');
% set(An_legend,'FontSize',30,'Textcolor','red');
% legend boxoff  
% hold on
fx=plot(xx,yy,'blue','LineWidth',3);
h_legend=legend(fx,'$f(x)$');
set(h_legend,'FontSize',30,'Textcolor','blue','interpreter','latex');
legend boxoff  
hold on
plot(x,y,'black*','LineWidth',4) 
hold on
plot(minpt_x,minpt_y,'red*','LineWidth',4)
hold on
z=x(1):1/1000:x(2);
plot(z,gn1(z),'green','LineWidth',2)
hold on
z=x(2):1/1000:x(3);
plot(z,gn2(z),'green','LineWidth',2)
hold on
z=x(3):1/1000:x(4);
plot(z,gn3(z),'green','LineWidth',2)
hold on
z=x(4):1/1000:x(5);
plot(z,gn4(z),'green','LineWidth',2)
hold on
z=x(5):1/1000:x(6);
plot(z,gn5(z),'green','LineWidth',2)
hold on
z=x(6):1/1000:x(7);
plot(z,gn6(z),'green','LineWidth',2)
axis([-3,3,-0.1,1.2])
axis off
