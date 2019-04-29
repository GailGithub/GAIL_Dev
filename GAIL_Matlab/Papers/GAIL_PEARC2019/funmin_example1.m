% funmin_g Example1

gail.InitializeDisplay
set(0,'defaultLineMarkerSize',20)

xAll = [];
fAll = [];
save fmin_ex1X xAll fAll

%Plot function
xplot = 0:0.001:1;
fplot = fmin_ex1(xplot);
h(1) = plot(xplot,fplot,'-');
set(h(1),'color',MATLABBlue)
hold on

%Plot fmin_g
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
[ffmg,outfmg] = funmin_g(@fmin_ex1,0,1);
h(2) = plot(mean(outfmg.intervals),ffmg,'.');
set(h(2),'color',MATLABGreen,'MarkerSize',80)
load fmin_ex1X xAll fAll
h(3) = plot(xAll,fAll,'.');
set(h(3),'color',MATLABGreen)

%Plot fminbnd
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
[xfmb,ffmb] = fminbnd(@fmin_ex1,0,1);
h(4) = plot(xfmb,ffmb,'.');
set(h(4),'color',MATLABOrange,'MarkerSize',80)
load fmin_ex1X xAll fAll
h(5) = plot(xAll,fAll,'.');
set(h(5),'color',MATLABOrange)

ylim([-6 1])
xlabel('$x$','interpreter','latex')
h_legend = legend([h(1) h(3) h(5)],{'$f(x)$','funmin\_g','fminbnd'},'Location','SouthEast');
set(h_legend,'interpreter','latex');
gail.save_eps('GAIL_PEARC2019_Output/','humps');



