% funmin_g Example1

gail.InitializeDisplay
set(0,'defaultLineMarkerSize',15)

xAll = [];
fAll = [];
save fmin_ex1X xAll fAll

%Plot function
xplot = 0:0.001:1;
fplot = fmin_ex1(xplot);
h(1) = plot(xplot,fplot,'-');
set(h(1),'color',MATLABBlue)
hold on

%Plot funmin_g
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
tic, [ffmg,outfmg] = funmin_g(@fmin_ex1,0,1); t1 = toc
h(2) = plot(mean(outfmg.intervals),ffmg,'.');
set(h(2),'color',MATLABGreen,'MarkerSize',80)
load fmin_ex1X xAll fAll
h(3) = plot(xAll,fAll,'.');
set(h(3),'color',MATLABGreen)

%Plot fminbnd
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
options = optimset('TolX',outfmg.abstol,'TolFun',outfmg.abstol); t2 = toc
tic, [xfmb,ffmb] = fminbnd(@fmin_ex1,0,1,options);
h(4) = plot(xfmb,ffmb,'.');
set(h(4),'color',MATLABOrange,'MarkerSize',80)
load fmin_ex1X xAll fAll
h(5) = plot(xAll,fAll,'.');
set(h(5),'color',MATLABOrange)

%Plot chebfun
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
tic, chebf = chebfun(@fmin_ex1,[0,1],'chebfuneps', outfmg.abstol, 'splitting','on'); t3=toc
chebfval = min(chebf); 
chebxvals = roots(diff(chebf));
[v,i] = min(abs(fmin_ex1(chebxvals)-chebfval));
chebxval = chebxvals(i);
chebn = length(chebf)
h(6) = plot(chebxval,chebfval,'o');
set(h(6),'color',MATLABPurple,'MarkerSize',20)
load fmin_ex1X xAll fAll
h(7) = plot(xAll,fAll,'.');
set(h(7),'color',MATLABPurple)


ylim([-6 1])
xlabel('$x$','interpreter','latex')
h_legend = legend([h(1) h(2) h(4) h(6) h(3) h(5) h(7)],{'$f(x)$','funmin\_g''s min',...
    'fminbnd''s min','chebfun''s min','funmin\_g''s sample','fminbnd''s sample','chebfun''s sample'},...
    'Location','Southeast');
set(h_legend,'interpreter','latex');

save fmin_ex1X ffmg outfmg xfmb ffmb chebxval chebfval xAll fAll
