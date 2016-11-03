%% Testing different routines for integration and function approximation

%% Initialization
format compact
close all %close all figures
clearvars %clear all variables
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels

%% Test Functions
% These are some test functions that we use to demonstrate the strengths
% and the weaknesses of the integration and function approximation
% routines.

%f=@(x) x; exactinteg=1/2;
%f=@(x) x.*x; exactinteg=1/3;
peaky=@(x,t,h) exp(-((x-t)/h).^2)/(0.5*h*sqrt(pi)*(erf((1-t)/h)+erf(t/h))); 
   %Gaussian peak
f = @(x) peaky(x,0.5,0.2); exactinteg = 1;

%% Plot the test function
xplot = 0:0.002:1;
plot(xplot,f(xplot),'-','linewidth',3)
xlabel('$x$')
ylabel('$f(x)$')


%% Integration routines
% Let's integrate the test function, i.e., we will compute
%
% \[ \int_0^1 f(x) \, {\rm d} x. \]
%
% The first method is MATLAB's built-in integration routine, which is based
% on L. F. Shampine, Vectorized Adaptive Quadrature in MATLAB, _Journal of
% Computational and Applied Mathematics_, *211*, 2008, pp. 131-140.

tic, integ1=integral (f,0,1); toc
err1 = abs(exactinteg-integ1)

%%
% The second method is the GAIL routine that has a theoretical
% justification

abstol = 1e-10;
tic, integ2=integral_g (f,0,1,abstol); toc 
err2 = abs(exactinteg-integ2) % gail routine

%%
% The third method is an improvement on the GAIL routine, which was
% proposed by

tic, integ3=integralNoPenalty_g (f,0,1,abstol); toc 
err3 = abs(exactinteg-integ3) % fred's routine

%%
% The fourth method is the Chebfun MATLAB toolbox

fcheb = chebfun(f,[0 1]); integ4=sum(fcheb); err4 = abs(exactinteg-integ4) % chebfun routine

%% Function Approx routines
xtest=sort(rand(1e6,1));
ftest=f(xtest);
f1=funappxglobal_g (f,0,1,abstol); errf1 = max(abs(ftest-ppval(f1,xtest))) % gail routine
