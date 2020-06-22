%% funmin_g
% 1-D guaranteed locally adaptive function optimization on [a,b]
%% Syntax
% fmin = *funmin_g*(f)
%
% fmin = *funmin_g*(f,a,b,abstol)
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol)
%
% fmin = *funmin_g*(f,in_param)
%
% [fmin, out_param] = *funmin_g*(f,...)
%% Description
%
% fmin = *funmin_g*(f) finds minimum value of function f on the
%  default interval [0,1] within the guaranteed absolute error tolerance
%  of 1e-6. Input f is a function handle.
%
% fmin = *funmin_g*(f,a,b,abstol) finds minimum value of
%  function f with ordered input parameters that define the finite
%  interval [a,b], and a guaranteed absolute error tolerance abstol.
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol) finds minimum
%  value of function f on the interval [a,b] with a guaranteed absolute
%  error tolerance. All three field-value pairs are optional and can be
%  supplied in different order.
%
% fmin = *funmin_g*(f,in_param) finds minimum value of function f
%  on the interval [in_param.a,in_param.b] with a guaranteed absolute
%  error tolerance in_param.abstol. If a field is not specified, the
%  default value is used.
%
% [fmin, out_param] = *funmin_g*(f,...) returns minimum value fmin
%  of function f and an output structure out_param.
%
% *Input Arguments*
%
% * f --- input function.
%
% * in_param.a --- left end point of interval, default value is 0.
%
% * in_param.b --- right end point of interval, default value is 1.
%
% * in_param.abstol --- guaranteed absolute error tolerance, default
%  value is 1e-6.
%
% *Optional Input Arguments*
%
% * in_param.ninit --- initial number of subintervals. Default to 20.
%
% * in_param.nmax --- cost budget, default value is 1e7.
%
% * in_param.maxiter --- max number of iterations, default value is 1000.
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
% * out_param.nmax --- cost budget
%
% * out_param.ninit --- initial number of points we use
%
% * out_param.npoints --- number of points needed to reach the guaranteed
%  absolute error tolerance
%
% <html>
% <ul type="square">
%  <li>out_param.exit --- this is a vector with two elements, for
%  tracking important warnings in the algorithm. The algorithm is
%  considered successful (with out_param.exit == [0 0]) if no flags arise
%  warning that the results are not guaranteed. The initial value is [0 0]
%  and the final value of this parameter is encoded as follows:</li>
%   <ul type="circle">
%    <li>[1 0]:  If reaching overbudget. It states whether
%                the max budget is attained without reaching the
%                guaranteed error tolerance.</li>
%    <li>[0 1]:  If reaching overiteration. It states whether
%                the max iterations is attained without reaching the
%                guaranteed error tolerance.</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.errest --- estimation of the absolute error bound
%
% * out_param.iter --- number of iterations
%
% * out_param.intervals --- the intervals containing point(s) where the
%  minimum occurs. Each column indicates one interval where the first raw
%  is the left point and the second row is the right point
%
%% Guarantee
%
% *Please check the details of the guarantee in [1].*
%
%% Examples
% *Example 1*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ with default input parameters.
  f = @(x) exp(0.01*(x-0.5).^2); [fmin,out_param] = funmin_g(f)


%%
% *Example 2*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ on $[-2,2]$ with error tolerance
% $10^{-7}$, cost budget $1000000$, initial number of points $10$.
f = @(x) exp(0.01*(x-0.5).^2);
[fmin,out_param] = funmin_g(f,-2,2,1e-7,10,1000000)



%%
% *Example 3*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ on $[-13,8]$ with error tolerance
% $10^{-7}$, cost budget $1000000$, initial number of points $100$.
clear in_param; in_param.a = -13; in_param.b = 8;
in_param.abstol = 1e-7;
in_param.ninit = 100;
in_param.nmax = 10^6;
[fmin,out_param] = funmin_g(f,in_param)


%%
% *Example 4*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ on $[-2,2]$ with error tolerance $10^{-5}$,
% cost budget $1000000$, initial number of points $64$.
f=@(x) exp(0.01*(x-0.5).^2);
[fmin,out_param] = funmin_g(f,'a',-2,'b',2,'ninit',64,'nmax',1e6,'abstol',1e-5)


%% See Also
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/fminbnd.html">fminbnd</a>
% </html>
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
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
% Adaption for Approximation and Minimization of Univariate Functions,"
% Journal of Complexity 40, pp. 17-33, 2017.
%
% [2] Xin Tong. A Guaranteed, "Adaptive, Automatic Algorithm for
% Univariate Function Minimization," MS thesis, Illinois Institute of
% Technology, 2014.
%
% [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [5] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
