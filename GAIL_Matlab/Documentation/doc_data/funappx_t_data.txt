%% Guarantee
%    
% |If the function to be approximated,|  $f$  |satisfies the cone condition|
%
% $$\|f''\|_\infty \le \frac { 2\mathrm{nstar} }{b-a } \left\|f'-\frac{f(b)-f(a)}{b-a}\right\|_\infty,$$
% 
% |then the|  $pp$  |output by this algorithm is guaranteed to
% satisfy|
%
% $$\| f-ppval(pp,)\|_{\infty} \le \mathrm{abstol},$$
%
% |provided the flag| $\mathrm{exceedbudget} = 0.$
%
% |And the upper bound of the cost is|
%
% $$\sqrt{ \frac{\mathrm{nstar}(b-a)^2 
% \|f''\|_\infty}{2 \times \mathrm{abstol}}} + 2 \times \mathrm{nstar}+4.$$
%
%
%% Examples
% *Example 1*

f = @(x) x.^2; [pp, out_param] = funappx_g(f)

% Approximate function x^2 with default input parameter to make the error
% less than 1e-6.
%%
% *Example 2*

[pp, out_param] = funappx_g(@(x) x.^2,0,100,1e-7,10,1000,1e8)

% Approximate function x^2 on [0,100] with error tolerence 1e-7, cost
% budget 10000000, lower bound of initial number of points 10 and upper
% bound of initial number of points 100

%%
% *Example 3*

clear in_param; in_param.a = -20; in_param.b = 20; in_param.nlo = 10;
in_param.nhi = 100; in_param.nmax = 1e8; in_param.abstol = 1e-7; 
[pp, out_param] = funappx_g(@(x) x.^2, in_param)

% Approximate function x^2 on [-20,20] with error tolerence 1e-7, cost
% budget 1000000, lower bound of initial number of points 10 and upper
% bound of initial number of points 100
%%
% *Example 4*

clear in_param; f = @(x) x.^2;
[pp, out_param] = funappx_g(f,'a',-10,'b',50,'nmax',1e6,'abstol',1e-8)

% Approximate function x^2 with error tolerence 1e-8, cost budget 1000000,
% lower bound of initial number of points 10 and upper
% bound of initial number of points 100
