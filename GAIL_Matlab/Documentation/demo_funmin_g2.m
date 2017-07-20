%% Compare *funmin_g* with *fminbnd*
% Author: Xin Tong, July 2017

%% Function definition
%
% Define a function with two minima as follows:
%
% \[ f(x) = -5 \exp(-100(x-0.2)^2) - \exp(-100(x-1)^2). \]
% 
close all; clearvars; format compact; format short;
f = @(x) -5*exp((-100*(x-0.2).^2))-exp((-100.*(x-1).^2));

%% Function minimization
% We use *funmin_g* to find the minimum of \(f\) over the interval
% \([a,b]\), where \(a = 0\) and \(b = 1.5\):
a = 0;
b = 1.5;
[fmin,out ] = funmin_g(f,a,b);
[xval,fval] = fminbnd(f,a,b);

%% Plot of the function and minima
% We plot \(f(x)\) and the global minimum value  returned
% by *funmin_g* and and a local minimum by *fminbnd* below:
figure;
x = a:1e-6:b;
fminvec = fmin.*ones(size(x));
plot(x,f(x),'r-',out.intervals,[fmin,fmin],'go',xval,fval,'b*'); 
ylim([-6 1])
xlabel('$x$','interpreter','latex')
h_legend=legend('$f(x)$','funmin\_g','fminbnd');
set(h_legend,'interpreter','latex');


%% References
%  
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%     _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%     Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%     GAIL: Guaranteed Automatic Integration Library (Version 2.2) [MATLAB
%     Software], 2017. Available from <http://gailgithub.github.io/GAIL_Dev/
%     GitHub>.

