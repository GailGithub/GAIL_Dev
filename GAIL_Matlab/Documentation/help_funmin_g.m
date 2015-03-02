%% funmin_g
% 1-D guaranteed global minimum value on [a,b] and the subset
% containing optimal solutions
%% Syntax
% fmin = *funmin_g*(f)
%
% fmin = *funmin_g*(f,a,b,abstol,TolX)
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol,'TolX',TolX)
%
% fmin = *funmin_g*(f,in_param)
%
% [fmin, out_param] = *funmin_g*(f,...)
%% Description
%
% fmin = *funmin_g*(f) finds minimum value of function f on the default
%  interval [0,1] within the guaranteed absolute error tolerance of 1e-6
%  and the X tolerance of 1e-3. Default initial number of points is 100
%  and default cost budget is 1e7. Input f is a function handle.
%
% fmin = *funmin_g*(f,a,b,abstol,TolX) finds minimum value of
%  function f with ordered input parameters that define the finite
%  interval [a,b], a guaranteed absolute error tolerance abstol and a
%  guaranteed X tolerance TolX.
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol,'TolX',TolX)
%  finds minimum value of function f on the interval [a,b] with a 
%  guaranteed absolute error tolerance abstol and a guaranteed X tolerance 
%  TolX. All five
%  field-value pairs are optional and can be supplied in different order.
%
% fmin = *funmin_g*(f,in_param) finds minimum value of function f on the
%  interval [in_param.a,in_param.b] with a guaranteed absolute error
%  tolerance in_param.abstol and a guranteed X tolerance in_param.TolX.
%  If a field is not specified, the default value is used.
%
% [fmin, out_param] = *funmin_g*(f,...) returns minimum value fmin of
%  function f and an output structure out_param.
%
% *Input Arguments*
%
% * f --- input function
%
% * in_param.a --- left end point of interval, default value is 0
%
% * in_param.b --- right end point of interval, default value is 1
%
% * in_param.abstol --- guaranteed absolute error tolerance, default
%  value is 1e-6.
%
% * in_param.TolX --- guaranteed X tolerance, default value is 1e-3.
%
% *Optional Input Arguments*
%
% * in_param.nlo --- lower bound of initial number of points we used,
%  default value is 10
%
% * in_param.nhi --- upper bound of initial number of points we used,
%  default value is 1000
%
% * in_param.nmax --- cost budget, default value is 1e7.
%
% *Output Arguments*
%
% * out_param.f --- input function
%
% * out_param.a --- left end point of interval
%
% * out_param.b --- right end point of interval
%
% * out_param.abstol --- guaranteed absolute error tolerance
%
% * out_param.TolX --- guaranteed X tolerance
%
% * out_param.nlo --- a lower bound of initial number of points we use
%
% * out_param.nhi --- an upper bound of initial number of points we use
%
% * out_param.nmax --- cost budget
%
% * out_param.ninit --- initial number of points we use
%
% * out_param.tau --- latest value of tau
%
% * out_param.npoints --- number of points needed to reach the guaranteed
%  absolute error tolerance or the guaranteed X tolerance
%
% * out_param.exitflag --- the state of program when exiting
%           0  Success
%           1  Number of points used is greater than out_param.nmax
%
% * out_param.errest --- estimation of the absolute error bound
%
% * out_param.volumeX --- the volume of intervals containing the point(s)
%  where the minimum occurs
%
% * out_param.tauchange --- it is 1 if out_param.tau changes, otherwise
%  it is 0
%
% * out_param.intervals --- the intervals containing point(s) where the
%  minimum occurs. Each column indicates one interval where the first
%  row is the left point and the second row is the right point.
%
%% Guarantee
%    
% If the function to be minimized, $f$ satisfies the cone condition
%
% $$\|f''\|_\infty \le  \frac {\tau}{b-a}\left\|f'-\frac{f(b)-f(a)}{b-a}
% \right\|_\infty,$$
%      
% then the $\mathrm{fmin}$ output by this algorithm is guaranteed to
% satisfy
%
% $$| \min f-\mathrm{fmin}| \le \mathrm{abstol},$$
%
% or
%
%     $$\mathrm{volumeX} \le \mathrm{TolX},$$
%
% provided the flag $\mathrm{exitflag} = 0.$
%
%
%% Examples
% *Example 1*

f=@(x) (x-0.3).^2+1; [fmin,out_param] = funmin_g(f)

% Minimize function (x-0.3)^2+1 with default input parameter.
%%
% *Example 2*

f=@(x) (x-0.3).^2+1;
[fmin,out_param] = funmin_g(f,-2,2,1e-7,1e-4,10,10,1000000)

% Minimize function (x-0.3)^2+1 on [-2,2] with error tolerence 1e-4, X
% tolerance 1e-2, cost budget 1000000, lower bound of initial number of
% points 10 and upper bound of initial number of points 10

%%
% *Example 3*

clear in_param; in_param.a = -13; in_param.b = 8;
in_param.abstol = 1e-7; in_param.TolX = 1e-4;
in_param.nlo = 10; in_param.nhi = 100;
in_param.nmax = 10^6;
[fmin,out_param] = funmin_g(f,in_param)

% Minimize function (x-0.3)^2+1 on [-13,8] with error tolerence 1e-7, X
% tolerance 1e-4, cost budget 1000000, lower bound of initial number of
% points 10 and upper bound of initial number of points 100

%%
% *Example 4*

f=@(x) (x-0.3).^2+1;
[fmin,out_param] = funmin_g(f,'a',-2,'b',2,'nhi',100,'nlo',10,...
    'nmax',1e6,'abstol',1e-4,'TolX',1e-2)

% Minimize function (x-0.3)^2+1 on [-2,2] with error tolerence 1e-4, X
% tolerance 1e-2, cost budget 1000000, lower bound of initial number of
% points 10 and upper bound of initial number of points 100
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
%% References
% [1]  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for
% Univariate Function Minimization. MS thesis, Illinois Institute of 
% Technology, 2014.
% 
% [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou, 
% "GAIL: Guaranteed Automatic Integration Library (Version 2.1)"
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software", Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, DOI:
% http://dx.doi.org/10.5334/jors.bb, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Sou-Cheng T. Choi, "Summary of the First Workshop On Sustainable
% Software for Science: Practice And Experiences (WSSSPE1)", Journal of
% Open Research Software, Volume 2, Number 1, e6, pp. 1-21, DOI:
% http://dx.doi.org/10.5334/jors.an, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
