%% CURVE 
% Approximate a highly fluctuating curve using funappx_g.

%% Function definition
%
% Define a highly fluctuating function as follows:
%
% \[ f(x) = x^2 sin (\frac{2 \pi}{ x^2} ). \] 
% 
close all; clear all; format compact; format short;
f = @(x) x.^2 .* sin((2*pi)./x.^2);

%% Function approximation
% We use funappx_g to approximate \(f\) over the interval
% \([a,b]\):
a = 0.1;
b = 2.5;
[q,out] = funappx_g(f, a, b);

%% Plot of the approximant 
% We plot out \(q(x)\) below:
figure;
x = a:1e-6:b;
plot(x,q(x)); 
xlabel('$x$','interpreter','latex')
ylabel('approximant')
axis tight

%% Plot of the apprroximation errors  
% The following plot shows that all pointwise absolute errors are less than
% the default tolerance of \(10^{-6}\).
figure;
plot(x,abs(f(x)-q(x))); 
xlabel('$x$','interpreter','latex')
ylabel('absolute error')
axis tight

%% A slightly different example
% If we changes \(a\) to a smaller number such as \(10^{-2}\), then even if
% we relax the tolerance to \(10^{-5}\), funappx_g may still return an
% approximant that fails to meet the tolerance. The reason is that \(f\) on
% \((a,b)\) is no longer in the cone of functions conducive for successful
% approximation.
a = 1e-2;
abstol = 1e-5;
[q2,out2] = funappx_g(f, a, b, abstol);
figure;
x = a:1e-6:b;
plot(x,abs(f(x)-q2(x))); 
xlabel('$x$','interpreter','latex')
ylabel('absolute error')
axis tight

%% A fix
% We can widen the cone by increasing the number of initial points given to
% funappx_g. 
inparam.a = a;
inparam.b = b;
inparam.abstol = abstol;
inparam.ninit = 5e6; 
inparam.nmax = inparam.ninit*10; 
[q3,out3] = funappx_g(f, inparam);
x = a:1.0/(out3.npoints*2):b;
figure;
plot(x,abs(f(x)-q3(x)));
xlabel('$x$','interpreter','latex')
ylabel('absolute error')
axis tight
 



 
