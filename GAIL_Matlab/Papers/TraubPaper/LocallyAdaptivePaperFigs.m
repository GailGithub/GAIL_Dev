function [funappx_g_success_rate, funmin_g_success_rate] = LocallyAdaptivePaperFigs(funappxRes,funminRes,colorfig)
% LOCALLYADPATIVEPAPERFIGS creates all figures and tables in the paper on local
% adpation for function approximation and optimization published in the
% Joseph Traub memorial issue in the Journal of Complexity:
%
%   Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%   Adaption for Approximation and Minimization of Univariate Functions,"
%   Journal of Complexity 40, pp. 17-33, 2017.
%
%
% Example 1: To plot all figures in the paper, run
%
%   LocallyAdaptivePaperFigs
%
% Example 2: To plot all figures and print summary tables from the input
% mat files (produced by running 'funappx_g_test' and 'funmin_g_test'), run
%
%   LocallyAdaptivePaperFigs('funappx_g_test-2016-10-17-18-20-32.mat', ...
%                            'funmin_g_test-2016-10-17-18-21-21.mat')

gail.InitializeDisplay %initialize the workspace and the display parameters
set(0,'defaultLineLineWidth',4) %thicker lines
if nargin < 3
   colorfig = true; %color figures
   if nargin < 2
      funminRes = []; %file for funmin_g testing
      if nargin < 1
         funappxRes = []; %file for funappx_g testing
      end
   end
end
if ~colorfig %black and white figures
   MATLABBlue = zeros(1,3);
   MATLABOrange = zeros(1,3);
   MATLABGreen = zeros(1,3);
end
whichdir = 'TraubPaperOutput';

%% Define function
f = @(x) x./(1+x.^3);
fpp = @(x) (6*x.^2 .*(-2 + x.^3))./(1 + x.^3).^3;
xnodes = 0:0.2:1;
fnodes = f(xnodes);
nnodes = numel(xnodes);
extra = 0.15;
xmin = min(xnodes) - extra;
xmax = max(xnodes) + extra;
xplot = xmin + (0:1e-3:1)*(xmax - xmin);
fppplot = fpp(xplot);
fppnodes = fpp(xnodes);
fppmax = max(abs(fppplot));
ylabelval = - 0.4;
yaxislvl = -0.2;

%% Plot second derivative of function
gail.RemovePlotAxes
plot([xnodes; xnodes],yaxislvl+0.1*[-ones(1,nnodes); ones(1,nnodes)],'-k')
plot([xmin xmax],yaxislvl*ones(1,2),'-k')
text(xnodes(1)-0.09,yaxislvl + ylabelval-0.02,'\(\beta - h_-\)');
text(xnodes(3)-0.03,yaxislvl + ylabelval-0.02,'\(\alpha\)');
text(xnodes(4)-0.03,yaxislvl + ylabelval-0.02,'\(\beta\)');
text(xnodes(6)-0.09,yaxislvl + ylabelval-0.02,'\(\alpha + h_+\)');
pbaspect([1 0.25 1])
fdd1 = diff(diff(fnodes(1:3))./diff(xnodes(1:3)))/diff(xnodes([1 3]));
fdd2 = diff(diff(fnodes(4:6))./diff(xnodes(4:6)))/diff(xnodes([4 6]));
plot(xplot,abs(fppplot),'-','color',MATLABBlue);
plot(xnodes, abs(fppnodes),'.','color',MATLABBlue)
plot(xnodes([1 3]),2*abs(fdd1)*ones(1,2), '-.', ...
   xnodes([4 6]),2*abs(fdd2)*ones(1,2), ...
   '-.', 'color',MATLABOrange);
plot(xnodes([1 3]),abs(fppnodes(1))*ones(1,2),':', ...
   xnodes([4 6]),abs(fppnodes(6))*ones(1,2), ...
   ':', 'color',MATLABPurple);
plot(xnodes([3 4]),abs(fppnodes(4))*ones(1,2), ...
   '--','color',MATLABGreen);
axis([xmin xmax yaxislvl-0.1 fppmax])
gail.save_eps(whichdir,'ConeDefineFig');

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
f3pp_xplot = f3pp(xplot);
[h,hLine1,hLine2] = plotyy(xplot,f3(xplot),xplot,f3pp_xplot);
set(hLine1, 'LineStyle','-','color',MATLABBlue)
set(hLine2, 'LineStyle','--','color',MATLABOrange)
set(h(1), 'YLim', [0, 1],'YTick',[0, 1])
set(h(2), 'YLim', [min(f3pp_xplot)-2, max(f3pp_xplot)+2],'YTick',[min(f3pp_xplot), max(f3pp_xplot)])
xlabel('\(x\)')
legend(h(1),{'\(f_1(x)\)', '\(f_1''''(x)\)'},'location', 'northeast','box','off')
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
in_param.ninit = 20;
in_param.abstol = 0.02;
in_param.a=-1;
in_param.b=1;
in_param.output_x = true;
c = -0.2;
delta = 0.3;
f3 = @(x) f3param(x,delta,c);
figure
[~,fappxout] = funappx_g(@(x) -f3(x),in_param);
newline = char(13);
disp([newline newline 'funappx_g used ' int2str(fappxout.npoints) ' points and ' ...
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
[~,fminout] = funmin_g(@(x) -f3(x), in_param);
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

%% Chebfun example
if exist('chebfun','file')
    delta = 0.2;
    f = @(x) f3param(x,delta,c);
    a = - 1;
    b = 1;
    chebfuntol = 1e-12;
    t0 = tic;
    chf = chebfun(f,[a,b],'chebfuneps', chebfuntol,'splitting','on');
    t_chebfun = toc(t0);
    x = a:0.00001:b;
    err = abs(chf(x) - f(x));
    [~,ind] = find(err > chebfuntol*10);
    figure;
    semilogy(x, err, '-','color',MATLABBlue);
    hold on
    semilogy(x(ind), err(ind), '.','color',MATLABOrange);
    small = max(-20,log10(0.1*min(err)));
    large = log10(10*max(err));
    axis([a b 10^small 10^large])
    xlabel('\(x\)')
    ylabel('Chebfun error')
    set(gca,'ytick',10.^(5*ceil(small/5):5:5*floor(large/5)))
    gail.save_eps('TraubPaperOutput', 'chebfun_errors');

    
    %% Time and error for funappx_g
    t0 = tic;
    [fhat, out] = funappx_g(f,a,b,chebfuntol);
    t_funappx = toc(t0);
    t_ratio = t_funappx/t_chebfun;
    err_funappx = max(abs(f(x) - fhat(x)));
end

%% Output funappx Table
funappx_g_success_rate = [];
if numel(funappxRes) > 0
   sorted_timeratio = 0;
   sorted_npointsratio = 0;
   trueerrormat = 0;
   exceedmat = 0;
   npoints = 0;
   time = 0;
   load(funappxRes)
   [fileID, fullPath] = gail.open_txt('TraubPaperOutput', 'funappx_g_test');
   fprintf(fileID,'\n');
   fprintf(fileID,'# of replications = %1.0f\n',nrep);
   fprintf(fileID,'   Test         Number of Points                    Time Used                          Success (%%)                                  Failure (%%)\n');
   fprintf(fileID,'  Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------\n');
   fprintf(fileID,'             Local      Global    Chebfun    Local       Global      Chebfun     Local        Global         Chebfun       Local         Global        Chebfun\n');
   fprintf(fileID,'                                                                                 No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn\n');
   npointslgratio = zeros(1,n);
   timelgratio = zeros(1,n);

   funmin_g_success_rate = zeros(length(permuted_index),1);   
   for i = permuted_index
     funappx_g_success_rate(i) = 100.0*sum(trueerrormat(i,1,:)<=abstol)/nrep;
      
     fprintf(fileID,'%9.0f %9.0f %9.0f  %9.0f %11.4f  %11.4f %11.4f  %6.0f %6.0f %6.0f %6.0f %6.0f   %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n',...
       [i mean(npoints(i,1,:)) mean(npoints(i,2,:)) mean(npoints(i,3,:))...
       mean(time(i,1,:)) mean(time(i,2,:)) mean(time(i,3,:))...
       100.0*sum(trueerrormat(i,1,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,1,:)<=abstol & (exceedmat(i,1,:)))/nrep ...
       100.0*sum(trueerrormat(i,2,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,2,:)<=abstol & (exceedmat(i,2,:)))/nrep ...
       100.0*sum(trueerrormat(i,3,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,3,:)<=abstol & (exceedmat(i,3,:)))/nrep...
       100.0*sum(trueerrormat(i,1,:)>abstol)/nrep  100.0*sum(trueerrormat(i,1,:)>abstol & (exceedmat(i,1,:)))/nrep ...
       100.0*sum(trueerrormat(i,2,:)>abstol)/nrep  100.0*sum(trueerrormat(i,2,:)>abstol & (exceedmat(i,2,:)))/nrep ...
       100.0*sum(trueerrormat(i,3,:)>abstol)/nrep  100.0*sum(trueerrormat(i,3,:)>abstol & (exceedmat(i,3,:)))/nrep]);
     npointslgratio(i) = mean(npoints(i,1,:))/mean(npoints(i,2,:));
     timelgratio(i) = mean(time(i,1,:))/mean(time(i,2,:));
   end
   fclose(fileID);
   type(fullPath)

%% Chebfun/funppx_g ratios
   figure
   t = ((1:nrep*n) -1/2)/(nrep*n);
   h = semilogx(sorted_timeratio(2,:),t,'color',MATLABBlue); hold on
   h = [h; semilogx(sorted_npointsratio(2,:),t,'--','color',MATLABOrange)]; 
       xlabel('Ratios'); 
       ylabel('Probability')
   legend(h,{'Time','\# Samples'},'location','northwest','box','off')
   hold off
   gail.save_eps('TraubPaperOutput', [algoname,'_test']);
end
 
%% Output funmin_g Table
if numel(funminRes) > 0
   trueerrormat = 0;
   exceedmat = 0;
   npoints = 0;
   time = 0;
   load(funminRes) 
   [fileID, fullPath] = gail.open_txt('TraubPaperOutput', 'funmin_g_test');
   fprintf(fileID,'\n');
   fprintf(fileID,'# of replications = %1.0f\n',nrep);
   fprintf(fileID,'   Test         Number of Points                    Time Used                          Success (%%)                                  Failure (%%)\n');
   fprintf(fileID,'  Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------\n');
   fprintf(fileID,'             funmin_g   fminbnd   Chebfun    funmin_g     fminbnd    Chebfun     funmin_g        fminbnd        Chebfun   funmin_g        fminbnd       Chebfun\n');
   fprintf(fileID,'                                                                                 No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn\n');
   npointslgratio = zeros(1,n);
   timelgratio = zeros(1,n);
   
   for i = permuted_index
       funmin_g_success_rate(i) = 100.0*sum(trueerrormat(i,1,:)<=abstol)/nrep;
       
       fprintf(fileID,'%9.0f %9.0f %9.0f  %9.0f %11.4f  %11.4f %11.4f  %6.0f %6.0f %6.0f %6.0f %6.0f   %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f \n',...
           [i mean(npoints(i,1,:)) mean(npoints(i,2,:)) mean(npoints(i,3,:))...
           mean(time(i,1,:)) mean(time(i,2,:)) mean(time(i,3,:))...
           100.0*sum(trueerrormat(i,1,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,1,:)<=abstol & (exceedmat(i,1,:)))/nrep ...
           100.0*sum(trueerrormat(i,2,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,2,:)<=abstol & (exceedmat(i,2,:)))/nrep ...
           100.0*sum(trueerrormat(i,3,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,3,:)<=abstol & (exceedmat(i,3,:)))/nrep...
           100.0*sum(trueerrormat(i,1,:)>abstol)/nrep  100.0*sum(trueerrormat(i,1,:)>abstol & (exceedmat(i,1,:)))/nrep ...
           100.0*sum(trueerrormat(i,2,:)>abstol)/nrep  100.0*sum(trueerrormat(i,2,:)>abstol & (exceedmat(i,2,:)))/nrep ...
           100.0*sum(trueerrormat(i,3,:)>abstol)/nrep  100.0*sum(trueerrormat(i,3,:)>abstol & (exceedmat(i,3,:)))/nrep]);
       npointslgratio(i) = mean(npoints(i,1,:))/mean(npoints(i,2,:));
       timelgratio(i) = mean(time(i,1,:))/mean(time(i,2,:));
   end
   fclose(fileID);
   type(fullPath)
end

%% SAMPLE OUTPUT

% >> LocallyAdaptivePaperFigs
% funappx_g used 65 points and 3 iterations
%    to satisfy an absolute error tolerance of 0.02
% the mesh width ranges from 0.025 to 0.1
% funmin_g used 43 points and 3 iterations
%    to satisfy an absolute error tolerance of 0.02
% the mesh width ranges from 0.025 to 0.1
% t_chebfun =
%     0.8860
% fhat = 
%   griddedInterpolant with properties:
% 
%             GridVectors: {[1x8388699 double]}
%                  Values: [1x8388699 double]
%                  Method: 'linear'
%     ExtrapolationMethod: 'linear'
% out = 
%            a: -1
%       abstol: 1.0000e-12
%            b: 1
%            f: @(x)f3param(x,delta,c)
%      maxiter: 1000
%        ninit: 20
%         nmax: 10000000
%     exitflag: [0 0 0 0 0]
%         iter: 21
%      npoints: 8388699
%       errest: 2.9057e-13
% t_funappx =
%     3.2602
% t_ratio =
%     3.6797
% err_funappx =
%    2.9775e-14

% >>  LocallyAdaptivePaperFigs('funappx_g_test-2016-10-17-18-20-32.mat', ...
%                              'funmin_g_test-2016-10-17-18-21-21.mat')
% funappx_g used 65 points and 3 iterations
%    to satisfy an absolute error tolerance of 0.02
% the mesh width ranges from 0.025 to 0.1
% funmin_g used 43 points and 3 iterations
%    to satisfy an absolute error tolerance of 0.02
% the mesh width ranges from 0.025 to 0.1
% t_chebfun =
%     0.1685
% fhat = 
%   griddedInterpolant with properties:
% 
%             GridVectors: {[1x8388699 double]}
%                  Values: [1x8388699 double]
%                  Method: 'linear'
%     ExtrapolationMethod: 'linear'
% out = 
%            a: -1
%       abstol: 1.0000e-12
%            b: 1
%            f: @(x)f3param(x,delta,c)
%      maxiter: 1000
%        ninit: 20
%         nmax: 10000000
%     exitflag: [0 0 0 0 0]
%         iter: 21
%      npoints: 8388699
%       errest: 2.9057e-13
% t_funappx =
%     3.4043
% t_ratio =
%    20.2017
% err_funappx =
%    2.9775e-14
% 
% 
% # of replications = 1000
%    Test         Number of Points                    Time Used                          Success (%)                                  Failure (%)
%   Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------
%              Local      Global    Chebfun    Local       Global      Chebfun     Local        Global         Chebfun       Local         Global        Chebfun
%                                                                                  No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn
%         3      6577     46439        116      0.0030       0.0024      0.0216     100      0    100      0      0        0      0      0      0      0    100      0
%         1      5017     26265         43      0.0032       0.0029      0.0063     100      0    100      0      3        0      0      0      0      0     97      0
%         2     15698     97106         22      0.0050       0.0079      0.0039     100      0    100      0      3        0      0      0      0      0     97      0
% 
% 
% # of replications = 1000
%    Test         Number of Points                    Time Used                          Success (%)                                  Failure (%)
%   Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------
%              funmin_g   fminbnd   Chebfun    funmin_g     fminbnd    Chebfun     funmin_g        fminbnd        Chebfun   funmin_g        fminbnd       Chebfun
%                                                                                  No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn
%         3       318         8        116      0.0025       0.0007      0.0272     100      0    100      0     14        0      0      0      0      0     86      0 
%         1       233        22         43      0.0024       0.0010      0.0070     100      0     27      0     60        0      0      0     73      0     40      0 
%         2       316         9         22      0.0026       0.0008      0.0039     100      0    100      0     35        0      0      0      0      0     65      0 