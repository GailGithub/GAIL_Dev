%% Guarantee
%    
% |If the function to be minimized,|  $f$  |satisfies the cone condition|
%
% $$\|f''\|_\infty \le  \tau\left\|f'-f(1)+f(0)\right\|_\infty,$$
% 
% |then the|  $\mathrm{fmin}$  |output by this algorithm is guaranteed to
% satisfy|
%
% $$| \min f-\mathrm{fmin}| \le \mathrm{abstol},$$
%
% |provided the flag| $\mathrm{exceedbudget} = 0.$
%
%
%% Examples
% *Example 1*

f = @(x) (x-0.3).^2+1; [fmin,out_param] = funmin_g(f)

% Minimize function (x-0.3)^2+1 with default input parameter.
%%
% *Example 2*

[fmin,out_param] = funmin_g(@(x) (x-0.3).^2+1,1e-4,1e-2,10,1e6)

% Minimize function (x-0.3)^2+1 with error tolerence 1e-4, X tolerance
% 1e-2, cost budget 1000000, initial number of points 10

%%
% *Example 3*

clear in_param; in_param.abstol = 1e-8; in_param.ninit = 10; 
in_param.nmax = 1e6; in_param.TolX = 1e-4;
[fmin,out_param] = funmin_g(@(x) (x-0.3).^2+1,in_param)

% Minimize function (x-0.3)^2+1 with error tolerence 1e-8, X tolerance
% 1e-4, cost budget 1000000, initial number of points 10

%%
% *Example 4*

clear in_param; f = @(x) (x-0.3)^2+1;
[fmin,out_param] = funmin_g(f,'ninit',52,'nmax',1e7,'abstol',1e-6,'TolX',1e-3)

% Minimize function (x-0.3)^2+1 with error tolerence 1e-6, X tolerance
% 1e-3, cost budget 10000000, initial number of points 52
