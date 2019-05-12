% min_pearc.m

if ~exist('chebfun','file') 
   warning('Chebfun is not installed.')
   return
end

%% Example 1a
close all; clearvars; format compact; format short;
funmin_example1;

load fmin_ex1X 
truefmin = fmin_ex1(0.15);
truexmin = 0.15;
funmin_g_max_abs_error = [max(abs(truexmin-mean(outfmg.intervals))), max(abs(truefmin-ffmg))]
fminbnd_max_abs_error = [max(abs(truexmin-xfmb)), max(abs(truefmin-ffmb))]
cheb_max_abs_error = [max(abs(truexmin-chebxval)), max(abs(truefmin-chebfval))]


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

tic, [fmin,out] = funmin_g(f, inparam); t1 = toc
options = optimset('TolX',out.abstol,'TolFun',out.abstol);
tic, [xval,fval,exitflag,output] = fminbnd(f,a,b, options); t2 = toc
tic, chebf = chebfun(f,[a,b],'chebfuneps', out.abstol, 'splitting','on'); t3 = toc
chebfval = min(chebf); 
chebxvals = roots(diff(chebf)); 
[v,i] = min(abs(f(chebxvals)-chebfval)); 
chebxval = chebxvals(i);
chebn = length(chebf)

truefmin = -2.99843616266006;
truexmin = 1.99843665971919;
funmin_g_max_abs_error = [max(abs(truexmin-mean(out.intervals))), max(abs(truefmin-fmin))]
fminbnd_max_abs_error = [max(abs(truexmin-xval)), max(abs(truefmin-fval))]
cheb_max_abs_error = [max(abs(truexmin-chebxval)), max(abs(truefmin-chebfval))]

x = 1.2:1e-6:b;
fminvec = fmin.*ones(size(x));
plot(x,f(x),'r-')
hold on
plot (out.intervals,[fmin,fmin],'go','MarkerSize',12)
hold on
plot(xval,fval,'b*','MarkerSize',12);
hold on
plot(chebxval,chebfval,'kd','MarkerSize',8);

ylim([-6 1])
xlabel('$x$','interpreter','latex')
h_legend = legend('$f(x)$','funmin\_g','fminbnd','chebfun','Location','SouthEast');
set(h_legend,'interpreter','latex');
gail.save_eps('GAIL_PEARC2019_Output/','sine');
