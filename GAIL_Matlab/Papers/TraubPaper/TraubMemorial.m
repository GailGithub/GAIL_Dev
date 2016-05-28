% function TraubMemorial
% TRAUBMEMORIAL create some figures for Traub memorial issue paper
% This script creates figures for 
%
% S.-C. T. Choi, Y. Ding, F. J. Hickernell, X. Tong, " ...

InitializeDisplay %initialize the workspace and the display parameters
MATLABBlue = [0, 0.447, 0.741];
MATLABOrange = [0.85,  0.325, 0.098];
MATLABPurple = [0.494,  0.184, 0.556];
MATLABGreen = [0.466,  0.674, 0.188];
%MATLABDkOrange = [0.85,  0.325, 0.098]*0.6;
%MATLABLtOrange = 0.5*[0.85,  0.325, 0.098] + 0.5*[1 1 1];
whichdir = 'TraubPaperOutput';

%% Sample functions with wildy oscillating second derivatives
f1 = @(x) x.^4 .* sin(1./((x==0)+x));
f1pp = @(x) (12*x.^2 - 1) .* sin(1./x) - 6*x .* cos(1./x);
f2 = @(x) f1(x) + 10.*x.^2;
f2pp = @(x) f1pp(x) + 20;
xplot = (-1:.001:1);
h = plot(xplot,f1(xplot),'-','color',MATLABBlue);
hold on
h=[h,plot(xplot,f2(xplot),'-','color',MATLABOrange)];
axis([-1 1 -1 11])
xlabel('\(x\)')
legend(h,{'\(f_1(x)\)', '\(f_2(x)\)'},'location', 'north','box','off')
gail.save_eps(whichdir, 'f1f2plot.eps');

figure
xplotclose = (-0.02:0.00001:0.02)';
h = plot(xplotclose,f1(xplotclose),'-','color',MATLABBlue);
%axis([-1 1 -1 11])
xlabel('\(x\)')
legend(h,{'\(f_1(x)\)'},'location', 'north','box','off')
print('-depsc',[whichdir 'f1closeplot.eps'])
gail.save_eps(whichdir, 'f1f2plot.eps');

figure
xplotclose = (-0.02:0.00001:0.02)';
h = plot(xplotclose,f2(xplotclose),'-','color',MATLABOrange);
axis([-0.02 0.02 -5e-4 4.5e-3])
xlabel('\(x\)')
legend(h,{'\(f_2(x)\)'},'location', 'north','box','off')
gail.save_eps(whichdir, 'f2closeplot.eps');

figure
h = plot(xplot,f1pp(xplot),'-',xplot,f2pp(xplot),'-');
axis([-1 1 -7 27])
xlabel('\(x\)')
legend(h,{'\(f''''_1(x)\)', '\(f''''_2(x)\)'},'location','best','box','off')
%print('-depsc',[whichdir 'f1ppf2ppplot.eps'])
gail.save_eps(whichdir, 'f1ppf2ppplot');

%% Sample functions with piecwise constant second derivatives
f3param = @(x,delta,c) (1./(2*delta.^2))*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
   - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta);
f3ppparam = @(x,delta,c) (1./(delta.^2))*(1 + sign(x-c-delta) - sign(x-c+delta)) ...
   .*(abs(x-c) <= 2*delta);
c = -0.2;
delta = 0.2;
f3 = @(x) f3param(x,delta,c);
f3pp = @(x) f3ppparam(x,delta,c);
figure
h = plot(xplot,f3(xplot),'-','color',MATLABPurple);
axis([-1 1 -0.2 1.2])
xlabel('\(x\)')
legend(h,{'\(f_3(x)\)'},'location', 'northeast','box','off')
gail.save_eps(whichdir, 'f3plot');

figure
h = plot(xplot,f3pp(xplot),'-','color',MATLABPurple);
%axis([-1 1 -0.2 1.2])
xlabel('\(x\)')
legend(h,{'\(f_3''''(x)\)'},'location', 'northeast','box','off')
print('-depsc',[whichdir 'f3ppplot.eps'])
gail.save_eps(whichdir, 'f3ppplot');

%% Lower Complexity Bounds Fooling Functions
rng(29);
n = 15;
%xpt = [-1 -1+2*sort(rand(1,n)) 1]; %points on [-1, 1]
xpt = [-1 (-1:2/n:1-2/n)+2*rand(1,n)/n 1]; %points on [-1, 1]
[~,wh] = max(diff(xpt));
c = (xpt(wh)+xpt(wh+1))/2;
delta = 1/(2*(n+1));
f3 = @(x) f3param(x,delta,c);
figure
h = plot(xplot,f3(xplot),'-','color',MATLABBlue);
hold on
h = [h; plot(xplot,-f3(xplot),'-','color',MATLABGreen)];
h = [h; plot(xpt(2:end-1),f3(xpt(2:end-1)),'.','color',MATLABOrange)];
axis([-1 1 -1.2 1.2])
xlabel('\(x\)')
legend(h,{'\(f_3(x)\)','\(-f_3(x)\)','\(\pm f_3(x_i)\)'}, ...
   'location', 'northeast','box','off')
gail.save_eps(whichdir, 'f3foolplot');

%% Sampling of hump functions for funappx and funmin
in_param.nlo = 10;
in_param.nhi = in_param.nlo;
in_param.abstol = 0.01;
in_param.a=-1;
in_param.b=1;
in_param.output_x = 1;
c = -0.2;
delta = 0.2;
f3 = @(x) f3param(x,delta,c);
figure
[~,fappxout] = funappx_g(@(x) -f3(x),in_param);
h = plot(xplot,-f3(xplot),'-', ...
   fappxout.x, fappxout.y,'.', ...
   fappxout.x, zeros(size(fappxout.x)),'.');
axis([-1 1 -1.2 0.2])
xlabel('\(x\)')
legend(h,{'\(-f_3(x)\)','\((x_i,-f_3(x_i))\)','\((x_i,0)\)'}, ...
   'location', 'east','box','off')
figure
[~,fappxout] = funmin_g(@(x) -f3(x),in_param);
h = plot(xplot,-f3(xplot),'-', ...
   fappxout.x, fappxout.y,'.', ...
   fappxout.x, zeros(size(fappxout.x)),'.');
axis([-1 1 -1.2 0.2])
xlabel('\(x\)')
legend(h,{'\(-f_3(x)\)','\((x_i,-f_3(x_i))\)','\((x_i,0)\)'}, ...
   'location', 'east','box','off')



