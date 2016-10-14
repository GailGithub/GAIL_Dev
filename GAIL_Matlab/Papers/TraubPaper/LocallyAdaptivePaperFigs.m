function LocallyAdaptivePaperFigs(colorfig)
% LOCALLYADPATIVEPAPERFIGS createS some figures for the paper on local
% adpation for function approximation and optimization submitted to the
% Joseph Traub memorial issue in the Journal of Complexity
%
% S.-C. T. Choi, Y. Ding, F. J. Hickernell, X. Tong, " ...

gail.InitializeDisplay %initialize the workspace and the display parameters
set(0,'defaultLineLineWidth',4) %thicker lines
whichdir = 'TraubPaperOutput';
if nargin < 1 
   colorfig = true; %color figures
end
if ~colorfig %black and white figures
   MATLABBlue = zeros(1,3);
   MATLABOrange = zeros(1,3);
   MATLABGreen = zeros(1,3);
end


%% Sample functions with piecwise constant second derivatives
xplot = (-1:.001:1);
f3param = @(x,delta,c) (1./(2*delta.^2))*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
   - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta);
f3ppparam = @(x,delta,c) (1./(delta.^2))*(1 + sign(x-c-delta) - sign(x-c+delta)) ...
   .*(abs(x-c) <= 2*delta);
c = -0.2;
delta = 0.2;
f3 = @(x) f3param(x,delta,c);
f3pp = @(x) f3ppparam(x,delta,c);

figure
[h,hLine1,hLine2] = plotyy(xplot,f3(xplot),xplot,f3pp(xplot));
set(hLine1, 'LineStyle','-','color',MATLABBlue)
set(hLine2, 'LineStyle','--','color',MATLABOrange)
set(h(1), 'YLim', [0, 1],'YTick',[0, 1])
set(h(2), 'YLim', [-30, 30],'YTick',[-30, 30])
xlabel('\(x\)')
legend(h(1),{'\(f_1(x)\)', '\(f_1''''(x)\)'},'location', 'northeast','box','off')
print('-depsc',[whichdir 'f3plot.eps'])
gail.save_eps(whichdir, 'f3plot');

%% Lower Complexity Bounds Fooling Functions
rng(29);
n = 15;
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
legend(h,{'\(f_1(x)\)','\(-f_1(x)\)','\(\pm f_1(x_i)\)'}, ...
   'location', 'northeast','box','off')
gail.save_eps(whichdir, 'f3foolplot');

%% Sampling of hump functions for funappx_g
in_param.nlo = 10;
in_param.nhi = in_param.nlo;
in_param.abstol = 0.02;
in_param.a=-1;
in_param.b=1;
in_param.output_x = 1;
c = -0.2;
delta = 0.3;
f3 = @(x) f3param(x,delta,c);
figure
[~,fappxout] = funappx_g(@(x) -f3(x),in_param);
disp(['funappx_g used ' int2str(fappxout.npoints) ' points and ' ...
   int2str(fappxout.iter) ' iterations'])
disp(['   to satisfy an absolute error tolerance of ' ...
   num2str(in_param.abstol)])
meshwidth = diff(fappxout.x);
disp(['the mesh width ranges from ' num2str(min(meshwidth)) ' to ' ...
   num2str(max(meshwidth))])
h = plot(xplot,-f3(xplot),'-','color',MATLABBlue);
hold on
h = [h; plot(fappxout.x, fappxout.y,'.','color',MATLABOrange)];
axis([-1 1 -1.2 0.2])
xlabel('\(x\)')
legend(h,{'\(-f_1(x)\)','\((x_i,-f_1(x_i))\)'}, ...
   'location', 'east','box','off')
gail.save_eps('TraubPaperOutput', 'sampling-funappxg');

%% Sampling of hump function for funmin_g
figure
[~,fminout] = funmin_g(@(x) -f3(x),in_param);
disp(['funmin_g used ' int2str(fminout.npoints) ' points and ' ...
   int2str(fminout.iter) ' iterations'])
disp(['   to satisfy an absolute error tolerance of ' ...
   num2str(in_param.abstol)])
meshwidth = diff(fminout.x);
disp(['the mesh width ranges from ' num2str(min(meshwidth)) ' to ' ...
   num2str(max(meshwidth))])
h = plot(xplot,-f3(xplot),'-','color',MATLABBlue);
hold on
h = [h; plot(fminout.x, fminout.y,'.','color',MATLABOrange)];
axis([-1 1 -1.2 0.2])
xlabel('\(x\)')
legend(h,{'\(-f_1(x)\)','\((x_i,-f_1(x_i))\)'}, ...
   'location', 'east','box','off')
gail.save_eps('TraubPaperOutput', 'sampling-funming');



