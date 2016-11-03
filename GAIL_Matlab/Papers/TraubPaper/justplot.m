% Just a sample not based on the real data

close all
%% initial curve
x=[-3,-2,-1,0,1,2,3];
y=[1.2,0.4,1,0,0.3,0.5,0.1];
xx = -3:1/1000:3;
yy = spline(x,y,xx);
yy1=interp1(x,y,xx);

%% create lower bound functions
minpt_x=[-3 -2 -2 -1 -1 0 0 1 1 2 2 3];
minpt_y=[0.6 -0.2 0 0.6 0.9 -0.1 -0.2 0.1 0.2 0.4 0.3 0];
gn1=@(x) -0.8.*x-1.8;
gn2=@(x) 0.6.*x+1.2;
gn3=@(x) -x-0.1; 
gn4=@(x) 0.3.*x-0.2; 
gn5=@(x) 0.2.*x; 
gn6=@(x) -0.3.*x+0.9; 


%% plot
hold on
plot(xx,zeros(1,size(xx,2)),'black','LineWidth',3); 
h_text=text(xx(end),0,'$U_n$');
set(h_text,'FontSize',30,'interpreter','latex');
hold on
An=plot(xx,yy1,'red');
An_legend=legend(An,'An(f)');
set(An_legend,'FontSize',30,'Textcolor','red');
legend boxoff  
hold on
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
axis([-3,3,-0.2,1.2])
axis off
