%% funappxab_g
% |One dimensional guaranteed function recovery on the interval [a,b].|
%% Syntax
% fappx = *funappxab_g*(f)
%
% fappx = *funappxab_g*(f,a,b,abstol,nlo,nhi,nmax)
%
% fappx = *funappxab_g*(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi,'nmax',nmax)
%
% fappx = *funappxab_g*(f,in_param)
%
% [fappx, out_param] = *funappxab_g*(f,...)
%% Description
% 
% fappx = *funappxab_g*(f) |recovers function|  f  |on the default interval
%   [a,b] by a piecewise linear interpolant fappx to within a guaranteed
%   absolute error of 1e-6. Default initial number of points is 52 and
%   default cost budget is 1e7.  Input| f |is a function handle. The
%   statement| y=f(x) |should accept a vector argument x and return a vector
%   y of function values that is the same size as x.|
%
% fappx = *funappxab_g*(f,a,b,abstol,nlo,nhi,nmax) |for given function|  f
%   |and the ordered input parameters with the interval a, b, guaranteed
%   absolute error abstol, lower bound of initial number of points nlo,
%   higher bound of initial number of points nhi and cost budget nmax. nlo
%   and nhi can be inputed as a vector or just one value as initial number
%   of points.|
%
% fappx =
% *funappxab_g*(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi,'nmax',nmax)
%   |recovers function|  f  |with the interval a, b, guaranteed absolute
%   error abstol, lower bound of initial number of points nlo, higher bound
%   of initial number of points nhi and cost budget nmax. All six
%   field-value pairs are optional and can be supplied in different order.|
%
% fappx = *funappxab_g*(f,in_param) |recovers function|  f  |with the
% interval in_param.a, in_param.b, guaranteed absolute error
% in_param.abstol, lower bound of initial number of points in_param.nlo,
% higher bound of initial number of points in_param.nhi and cost budget
% in_param.nmax. If a field is not specified, the default value is used.|
%
% [fappx, out_param] = *funappxab_g*(f,...) |returns function approximation
% fappx and an output structure out_param.|
%
% *Input Arguments*
% 
% * f --- |function handle|
%
% * in_param.a --- |left end point of interval, default value is 0|
%
% * in_param.b --- |right end point of interval, default value is 1|
%
% * in_param.abstol --- |guaranteed absolute error, default value
%  is 1e-6|
%
% * in_param.nlo --- |lower bound of initial number of points we used,
%  default value is 52|
%
% * in_param.nhi --- |higher bound of initial number of points we used,
%  default value is 52|
%
% * in_param.nmax --- |cost budget, default value is 1e7|
%
% *Output Arguments*
%
% * out_param.nmax --- |cost budget|
% 
% * out_param.exceedbudget  --- |it is 0 if the number of points used in the 
%   construction of fappx is less than cost budget, 1 otherwise.|
% 
% * out_param.ninit --- |initial number of points we used|
%
% * out_param.npoints --- |number of points we need to reach the guaranteed
% absolute error|
% 
% * out_param.errorbound --- |estimation of the approximation absolute error
% bound|
% 
% * out_param.tau --- |latest value of tau|
%
% * out_param.a --- |left end point of interval|
%
% * out_param.b --- |right end point of interval|
%
% * out_param.abstol --- |guaranteed absolute error|
% 
% * out_param.nlo --- |lower bound of initial number of points we used,
% default value is 52|
%
% * out_param.nhi --- |higher bound of initial number of points we used,
% default value is 52|
%
%% Guarantee
%    
% |If the function to be approximated,|  $f$  |satisfies the cone condition|
%
% $$\|f''\|_\infty \le \frac { \mathrm{out\_param.tau} }{b-a} \left\|f'-\frac{f(b)-f(a)}{b-a}\right\|_\infty,$$
% 
% |then the|  $fappx$  |output by this algorithm is guaranteed to
% satisfy|
%
% $$\| f-fappx \|_{\infty} \le \mathrm{out\_param.abstol}$$
%
% |provided the flag| 
%
% $$\mathrm{out\_param.exceedbudget} = 0.$$
%
%% Examples
% *Example 1*

f = @(x) x.^2; [fappx, out_param] = funappxab_g(f)

% Approximate function x^2 with default input parameter to make the error
% less than 1e-6.
%%
% *Example 2*

[fappx, out_param] = funappxab_g(@(x) x.^2,0,100,1e-7,52,52,1e8)

% Approximate function x^2 on [0,100] with error tolerence 1e-7, cost
% budget 10000000 and initial number of points 52

%%
% *Example 3*

clear in_param; in_param.a = -20; in_param.b = 20; in_param.nlo = 10;
in_param.nhi = 100; in_param.nmax = 1e8; in_param.abstol = 1e-7; 
[fappx, out_param] = funappxab_g(@(x) x.^2, in_param)

% Approximate function x^2 on [-20,20] with error tolerence 1e-7, cost
% budget 1000000, lower bound of initial number of points 10 and upper
% bound of initial number of points 100
%%
% *Example 4*

clear in_param; f = @(x) x.^2;
[fappx, out_param] = funappxab_g(f,'a',-10,'b',50,'nmax',1e6,'abstol',1e-8)

% Approximate function x^2 with error tolerence 1e-8, cost budget 1000000
% and initial number of points 10
%% See Also
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
%% Reference
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, The
% Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, Not Balls,
% Journal of Complexity 30 (2014) 21–45
%