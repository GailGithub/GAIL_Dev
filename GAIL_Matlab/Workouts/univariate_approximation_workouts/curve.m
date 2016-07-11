%% CURVE 
% Approximate a curve using funappxNoPenalty_g.

%% Function definition
%
% Define a highly fluctuating function as follows:
%
% \[ f(x) = x^2 sin(\frac{2 \pi}{ x^2} ). \] 
% 
f = @(x) x.^2 .* sin((2*pi)./x.^2);

%% Function approximation
% We use funappxNoPenalty_g to approximate \(f\) over the interval
% \([a,b]\):
a = 1e-2;
b = 2.5;
[q,out] = funappxNoPenalty_g(f, a, b);

%% Plot of the approximant 
% We plot out \(q(x)\) below:
figure;
x = a:0.0001:b;
plot(x,q(x)); 
axis tight

%% Plot of the apprroximation errors  
% The following plot shows that all pointwise absolute errors are less than the
% default tolerance of \(10^{-6}\).
figure;
plot(x,abs(f(x)-q(x))); 
axis tight
