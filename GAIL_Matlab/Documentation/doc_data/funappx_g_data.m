%% Guarantee
%
% For $[a,b]$, there exists a partition
%
% $$ P=\{[t_0,t_1], [t_1,t_2],  \ldots, [t_{L-1},t_L]\},  a=t_0 < t_1 < \cdots < t_L=b.$$
% 
% If the function to be approximated,  $f$ satisfies the cone condition
%
% $$\|f''\|_\infty \le \frac { 2\mathrm{nstar} }{t_l-t_{l-1} } \left\|f'-\frac{f(t_l)-f(t_{l-1})}{t_l-t_{l-1}}\right\|_\infty$$
% 
% for each sub interval $[t_{l-1},t_l]$, where $1 \le l \le L$,
% then the $pp$ |output by this algorithm is guaranteed to satisfy
%
% $$\| f-ppval(pp,)\|_{\infty} \le \mathrm{abstol}.$$
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
