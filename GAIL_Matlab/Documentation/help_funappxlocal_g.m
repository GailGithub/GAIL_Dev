%% funappxlocal_g
% |1-D guaranteed function recovery on closed interval [a,b].|
%% Syntax
% pp = *funappxlocal_g*(f)
%
% pp = *funappxlocal_g*(f,a,b,abstol,taulo,tauhi)
%
% pp = *funappxlocal_g*(f,'a',a,'b',b,'abstol',abstol,'taulo',taulo,'tauhi',tauhi)
%
% pp = *funappxlocal_g*(f,in_param)
%
% [pp, out_param] = *funappxlocal_g*(f,...)
%% Description
% 
% pp = *funappxlocal_g*(f) |recovers function|  f  |on the default interval
%  [0,1] by a piecewise polynomial structure pp to within the guaranteed
%  absolute error tolerance of 1e-6. Input| f |is a function handle. The
%  statement| y=f(x) |should accept a vector argument x and return a
%  vector y of function values that is the same size as x. Output pp
%   may be evaluated via| ppval.
%
% pp = *funappxlocal_g*(f,a,b,abstol,taulo,tauhi,nmax) |for a given
% function| f |and the ordered input parameters that define the finite 
% interval [a,b], a guaranteed absolute error tolerance bstol, a lower
% bound of cone condition taulo, and an upper bound of cone condition tauhi.|
%
% pp =
% *funappxlocal_g*(f,'a',a,'b',b,'abstol',abstol,'taulo',taulo,'tauhi',tauhi,'nmax',nmax)
%  |recovers function|  f  |on the finite interval [a, b], given a
%  guaranteed absolute error tolerance abstol,a lower bound of initial cone
%  condition taulo, and an upper bound of initial cone condition tauhi.
%  All five field-value pairs are optional and can be supplied in
%  different order.|
%
% pp = *funappxlocal_g*(f,in_param) |recovers function|  f  |on the finite
%  interval [in_param.a, in_param.b], given a guaranteed absolute error
%  tolerance in_param.abstol, a lower bound of initial cone condition
%  in_param.taulo, an upper bound of initial cone condition
%  in_param.tauhi. If a field is not specified, the default value is
%  used.|
%
% [pp, out_param] = *funappxlocal_g*(f,...) |returns a piecewise polynomial
%   structure pp and an output structure out_param.|
%
% *Input Arguments*
% 
% * f --- |function handle|
%
% * in_param.a --- |left end point of interval, default value is 0|
%
% * in_param.b --- |right end point of interval, default value is 1|
%
% * in_param.abstol --- |guaranteed absolute error tolerance, default value
%  is 1e-6|
%
% * in_param.taulo --- |lower bound of initial number of points we used,
%  default value is 9|
%
% * in_param.tauhi --- |upper bound of initial number of points we used,
%  default value is 100|
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
% * out_param.ninit --- |initial number of points we used|
%
% * out_param.npoints --- |number of points we need to reach the guaranteed
% absolute error|
% 
% * out_param.errorbound --- |an upper bound of the absolute error|
% 
% * out_param.tau --- |a vector indicate the cone condition of each
% subinterval|
%
% * out_param.a --- |left end point of interval|
%
% * out_param.b --- |right end point of interval|
%
% * out_param.abstol --- |guaranteed absolute error|
% 
% * out_param.taulo --- |a lower bound of cone condtion|
%
% * out_param.tauhi --- |an upper bound of cone condtion|
%
%% Guarantee

%% Examples
% *Example 1*

f = @(x) exp(-100*x.^2); [pp, out_param] = funappxlocal_g(f)

% Approximate function e^{-100x^2} with default input parameter to make the
% error less than 1e-6.
%%
% *Example 2*

[pp, out_param] = funappxlocal_g(@(x) exp(-100*x.^2),0,10,1e-7,10,1000)

% Approximate function e^{-100x^2} on [0,10] with error tolerence
% 1e-7,lower bound of initial number of points 10 and upper bound of
% initial number of points 100

%%
% *Example 3*

clear in_param; in_param.a = -10; in_param.b = 10; in_param.taulo = 10;
in_param.tauhi = 100; in_param.abstol = 1e-7; 
[pp, out_param] = funappxlocal_g(@(x) exp(-100*x.^2), in_param)

% Approximate function e^{-100x^2} on [-20,20] with error tolerence 1e-7,
% cost budget 1000000, lower bound of initial number of points 10 and upper
% bound of initial number of points 100
%% See Also
%
% <html> <a href="help_integral_g.html">integral_g</a> </html>
%
% <html> <a href="help_meanMC_g.html">meanMC_g</a> </html>
%
% <html> <a href="help_cubMC_g.html">cubMC_g</a> </html>
%
%% Reference
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, The
% Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, Not Balls,
% Journal of Complexity 30 (2014) 21–45
%
% [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
% Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
% 1.3.0)" [MATLAB Software], 2014. Available from
% http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
%