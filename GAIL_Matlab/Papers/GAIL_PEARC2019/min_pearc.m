% min_pearc.m


%% Example 1a
close all; clearvars; format compact; format short;
gail.InitializeDisplay %initialize the workspace and the display parameters
set(0,'defaultLineLineWidth',4) %thicker lines
f = @(x) -5*exp((-100*(x-0.2).^2))-exp((-100.*(x-1).^2));
a = 0;
b = 1.5;

tic, [fmin,out ] = funmin_g(f,a,b), t1 = toc
tic, [xval,fval,exitflag,output] = fminbnd(f,a,b), t2 = toc

truefmin = -5;
truexmin = 0.2;
funmin_g_max_abs_error = [max(abs(truexmin-mean(out.intervals))), max(abs(truefmin-fmin))]
fminbnd_max_abs_error = [max(abs(truexmin-xval)), max(abs(truefmin-fval))]

figure;
x = a:1e-6:b;
fminvec = fmin.*ones(size(x));

plot(x,f(x),'r-',out.intervals,[fmin,fmin],'go',xval,fval,'b*');
ylim([-6 1])
xlabel('$x$','interpreter','latex')
h_legend = legend('$f(x)$','funmin\_g','fminbnd','Location','SouthEast');
set(h_legend,'interpreter','latex');
gail.save_eps('GAIL_PEARC2019_Output/','humps');

%% Example 1b
figure; clearvars; format compact; format short;
gail.InitializeDisplay %initialize the workspace and the display parameters
set(0,'defaultLineLineWidth',4) %thicker lines
f = @(x) sin(10*pi*x.^4)-x;
a = 0; b = 2;

inparam.a = a;
inparam.b = b;
inparam.ninit = 1000;
inparam.nmax = inparam.ninit*10;

tic, [fmin,out] = funmin_g(f, inparam), t1 = toc
options = optimset('TolX',1e-6,'TolFun', 1e-6);
tic, [xval,fval,exitflag,output] = fminbnd(f,a,b, options), t2 = toc

truefmin = -2.99843616266006;
truexmin = 1.99843665971919;
funmin_g_max_abs_error = [max(abs(truexmin-mean(out.intervals))), max(abs(truefmin-fmin))]
fminbnd_max_abs_error = [max(abs(truexmin-xval)), max(abs(truefmin-fval))]

x = 1.2:1e-6:b;
fminvec = fmin.*ones(size(x));
plot(x,f(x),'r-',out.intervals,[fmin,fmin],'go',xval,fval,'b*');
ylim([-6 1])
xlabel('$x$','interpreter','latex')
h_legend = legend('$f(x)$','funmin\_g','fminbnd','Location','SouthEast');
set(h_legend,'interpreter','latex');
gail.save_eps('GAIL_PEARC2019_Output/','sine');
