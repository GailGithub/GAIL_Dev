%% funappx_g
% |One dimensional guaranteed function recovery on the interval [0,1].|
%% Syntax
% fappx = *funappx_g*(f)
%
% fappx = *funappx_g*(f,abstol,ninit,nmax)
%
% fappx = *funappx_g*(f,'abstol',abstol,'ninit',ninit,'nmax',nmax)
%
% fappx = *funappx_g*(f,in_param)
%
% [fappx, out_param] = *funappx_g*(f,...)
%% Description
% 
% fappx = *funappx_g*(f) |recovers function f on the interval [0,1] by a
% piecewise linear interpolant fappx to within a guaranteed absolute 
% error of 1e-6. Default initial number of points is 52 and default cost
% budget is 1e7. Input f is a function handle. The statement Y=f(X) should
% accept a vector argument X and return a vector Y of function values that
% is the same size as X.|
%
% fappx = *funappx_g*(f,abstol,ninit,nmax) |for given function f and the 
%   ordered input parameters with the guaranteed absolute error _abstol_,
%   initial number of points _ninit_ and cost budget _nmax_.|
%
% fappx = *funappx_g*(f,'abstol',abstol,'ninit',ninit,'nmax',nmax) |recovers 
%   function f with the guaranteed absolute error _abstol_, initial number
%   of points _ninit_, and cost budget _nmax_ . All three field-value pairs are
%   optional and can be supplied in different order.|
%
% fappx = *funappx_g*(f,in_param) |recovers function f with the guaranteed
% absolute error _in_param.abstol_, initial number of points _in_param.nint_,
% and cost budget _in_param.nmax_. If a field is not specified, the default
% value is used.|
%
% [fappx, out_param] = *funappx_g*(f,...) |returns an approximated function
% fappx and an output structure out_param.|
%
% *Input Arguments*
% 
% * f --- |function handle|
%
% * in_param.abstol --- |guaranteed absolute error, default value
%  is 1e-6|
%
% * in_param.ninit --- |initial number of points, default value is 52|
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
% * out_param.abstol --- |guaranteed absolute error|
% 

%% Examples
% *Example 1*

format short; f = @(x) x.^2; [fappx, out_param] = funappx_g(f)

% Approximate function x^2 with default input parameter to make the error
% less than 1e-6.
%%
% *Example 2*

format short; 
[fappx, out_param] = funappx_g(@(x) x.^2,1e-8,10,100000)

% Approximate function x^2 with error tolerence 1e-8, cost budget 100000
% and initial number of points 10

%%
% *Example 3*

format short; clear in_param; in_param.ninit = 10; in_param.nmax = 1e6; 
in_param.abstol = 1e-7; [fappx, out_param] = funappx_g(@(x) x.^2, in_param)

% Approximate function x^2 with error tolerence 1e-7, cost budget 1000000
% and initial number of points 10
%%
% *Example 4*

format short; clear in_param;
[fappx, out_param] = funappx_g(@(x) x.^2,'ninit',10,'nmax',1e6,'abstol',1e-8)

% Approximate function x^2 with error tolerence 1e-8, cost budget 1000000
% and initial number of points 10
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
%% Reference
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
% The Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, 
% Not Balls, submitted for publication, arXiv.org:1303.2412 
% [math.NA]}, 2013.
%