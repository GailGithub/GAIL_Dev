%% funappx_g
% |1-D guaranteed function recovery on a closed interval [a,b]|
%% Syntax
% pp = *funappx_g*(f)
%
% pp = *funappx_g*(f,a,b,abstol,nlo,nhi,nmax,maxiter)
%
% pp = *funappx_g*(f,in_param)
%
% [pp, out_param] = *funappx_g*(f,...)
%% Description
%
% pp = *funappx_g*(f) approximates function f on the default interval
%  [0,1] by a piecewise polynomial structure pp within the guaranteed
%  absolute error tolerance of 1e-6. Input f is a function handle. The
%  statement y = f(x) should accept a vector argument x and return a
%  vector y of function values that is of the same size as x. Output pp
%  may be evaluated via PPVAL.
%
% pp = *funappx_g*(f,a,b,abstol,nlo,nhi,nmax,maxiter) for a given
%  function f and the ordered input parameters that define the finite
%  interval [a,b], a guaranteed absolute error tolerance abstol, a lower
%  bound of initial number of points nlo, an upper bound of initial number
%  of points nhi, a cost budget nmax and max number of iteration maxiter.
%
%  pp =
%  FUNAPPX_G(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi,'nmax',nmax,'maxiter',maxiter)
%  recovers function f on the finite interval [a,b], given a guaranteed
%  absolute error tolerance abstol, a lower bound of initial number of
%  points nlo, an upper bound of initial number of points nhi, a cost
%  budget nmax and max number of iteration maxiter. All seven field-value
%  pairs are optional and can be supplied in different order.
%
% pp = *funappx_g*(f,in_param) recovers function f on the finite
%  interval [in_param.a,in_param.b], given a guaranteed absolute error
%  tolerance in_param.abstol, a lower bound of initial number of points
%  in_param.nlo, an upper bound of initial number of points in_param.nhi,
%  a cost budget in_param.nmax and max number of iteration
%  in_param.maxiter. If a field is not specified, the default value is
%  used.
%
% [pp, out_param] = *funappx_g*(f,...) returns a piecewise polynomial
%  structure pp and an output structure out_param.
%
% *Input Arguments*
%
% * in_param.a --- |left end point of interval, default value is 0|
%
% * in_param.b --- |right end point of interval, default value is 1|
%
% * in_param.abstol --- |guaranteed absolute error tolerance, default
%  value is 1e-6|
%
% * in_param.nlo --- |lower bound of initial number of points we used,
%  default value is 10|
%
% * in_param.nhi --- |upper bound of initial number of points we used,
%  default value is 1000|
%
% * in_param.nmax --- |when number of points hits the value, iteration
%  will stop, default value is 1e7|
%
% * in_param.maxiter --- |max number of interation, default value is 1000|
%
% *Output Arguments*
%
% * pp.form --- |pp means piecewise polynomials|
%
% * pp.breaks --- |show the location of interpolation points|
%
% * pp.coefs --- |coefficients for piecewise linear polynomials|
%
% * pp.pieces --- |number of piecewise linear polynomials|
%
% * pp.order --- |be 2 as we use piecewise linear polynomials|
%
% * pp.dim --- |be 1 as we do univariate approximation|
%
% * pp.orient --- |always be 'first'|
%
% * out_param.ninit --- |initial number of points we use for each sub
%  interval|
%
% * out_param.npoints --- |number of points we need to reach the
%  guaranteed absolute error tolerance|
%
% * out_param.errorbound --- |an upper bound of the absolute error|
%
% * out_param.nstar --- |final value of the parameter defining the cone of
%  functions for which this algorithm is guaranteed for each
%  subinterval; nstar = ninit-2 initially|
%
% * out_param.a --- |left end point of interval|
%
% * out_param.b --- |right end point of interval|
%
% * out_param.abstol --- |guaranteed absolute error tolerance|
%
% * out_param.nlo --- |a lower bound of initial number of points we use|
%
% * out_param.nhi --- |an upper bound of initial number of points we use|
%
% * out_param.nmax --- |when number of points hits the value, iteration
%  will stop|
%
% * out_param.maxiter --- |max number of interation|
%
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
%% See Also
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_funmin_g.html">funmin_g</a>
% </html>
%
%% References
%
% [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls, Journal of Complexity 30 (2014), pp. 21-45.
%
%
% [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% "GAIL: Guaranteed Automatic Integration Library (Version 2.0)" [MATLAB
% Software], 2014. Available from http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
%
