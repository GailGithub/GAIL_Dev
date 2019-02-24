%% Parameter settings
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',1) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',7) %latex axis labels

c1 = 0.1; c2 = 0.1;
alpha = 10;
beta = 5;
% f= @(X,Y) sqrt(X.*Y.^2)./(sin(1./(X+.05))+2);
f= @(X,Y) sin(((alpha*(X-c1)).^2+(beta*(Y-c2)).^2).^(0.6))./((alpha*(X-c1)).^2+(beta*(Y-c2)).^2).^(.3);
phi = @(x) 1-abs(2*x-1);
[X,Y] = meshgrid(0:.001:1);

mesh(X,Y,f(X,Y))
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
% title('$f(x_1,x_2)=\frac{\sin\left((100(x_1-0.1)^2+25(x_2-0.1)^2)^{0.6}\right)}{\left(100(x_1-0.1)^2+25(x_2-0.1)^2\right)^{0.3}}$','interpreter','latex')
title('$f(x_1,x_2)$','interpreter','latex')
axis([0 1 0 1 -1 1])
print('-depsc', strcat('baker_not.eps'))
figure
mesh(X,Y,f(phi(X),phi(Y)))
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
title('$f(\phi(x_1),\phi(x_2))$','interpreter','latex')
axis([0 1 0 1 -1 1])
print('-depsc', strcat('baker_with.eps'))