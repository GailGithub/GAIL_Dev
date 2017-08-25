%% integral_g
% 1-D guaranteed function integration using trapezoidal rule
%% Syntax
% q = *integral_g*(f)
%
% q = *integral_g*(f,a,b,abstol)
%
% q = *integral_g*(f,'a',a,'b',b,'abstol',abstol)
%
% q = *integral_g*(f,in_param)
%
% [q, out_param] = *integral_g*(f,...)
%% Description
%
% q = *integral_g*(f) computes q, the definite integral of function f on
%  the interval [a,b] by trapezoidal rule with in a guaranteed absolute
%  error of 1e-6. Default starting number of sample points taken is 100
%  and default cost budget is 1e7. Input f is a function handle. The
%  function y = f(x) should accept a vector argument x and return a vector
%  result y, the integrand evaluated at each element of x.
%
% q = *integral_g*(f,a,b,abstol) computes q, the definite integral of
%  function f on the finite interval [a,b] by trapezoidal rule with the
%  ordered input parameters, and guaranteed absolute error tolerance
%  abstol.
%
% q = *integral_g*(f,'a',a,'b',b,'abstol',abstol) computes q, the definite
%  integral of function f on the finite interval [a,b] by trapezoidal rule
%  within a guaranteed absolute error tolerance abstol. All four
%  field-value pairs are optional and can be supplied.
%
% q = *integral_g*(f,in_param) computes q, the definite integral of
%  function f by trapezoidal rule within a guaranteed absolute error
%  in_param.abstol. If a field is not specified, the default value is
%  used.
%
% [q, out_param] = *integral_g*(f,...) returns the approximated integration
%  q and output structure out_param.
%
% *Input Arguments*
%
% * f --- input function
%
% * in_param.a --- left end of the integral, default value is 0
%
% * in_param.b --- right end of the integral, default value is 1
%
% * in_param.abstol --- guaranteed absolute error tolerance, default value
%  is 1e-6
% 
% *Optional Input Arguments*
%
% * in_param.nlo --- lowest initial number of function values used, default
%  value is 10
%
% * in_param.nhi --- highest initial number of function values used,
%  default value is 1000
%
% * in_param.nmax --- cost budget (maximum number of function values),
%  default value is 1e7
%
% * in_param.maxiter --- max number of iterations, default value is 1000
% 
% *Output Arguments*
%
% * q --- approximated integral
%
% * out_param.f --- input function
%
% * out_param.a --- low end of the integral
%
% * out_param.b --- high end of the integral
%
% * out_param.abstol --- guaranteed absolute error tolerance
% 
% * out_param.nlo --- lowest initial number of function values
%
% * out_param.nhi --- highest initial number of function values
%
% * out_param.nmax --- cost budget (maximum number of function values)
%
% * out_param.maxiter --- max number of iterations
%
% * out_param.ninit --- initial number of points we use, computed by nlo
%  and nhi
%
% * out_param.exceedbudget --- it is true if the algorithm tries to use
%  more points than cost budget, false otherwise.
% 
% * out_param.tauchange --- it is true if the cone constant has been
%  changed, false otherwise. See [1] for details. If true, you may wish to
%  change the input in_param.ninit to a larger number.
% 
% * out_param.iter --- number of iterations
%
% * out_param.npoints --- number of points we need to 
%  reach the guaranteed absolute error tolerance abstol.
%
% * out_param.errest --- approximation error defined as the differences
%  between the true value and the approximated value of the integral.
%
% * out_param.nstar --- final value of the parameter defining the cone of
%  functions for which this algorithm is guaranteed; nstar = ninit-2
%  initially and is increased as necessary
%
%% Guarantee
%    
% If the function to be integrated, \(f\) satisfies the cone condition
%
% \[\|f''\|_1 \le \frac { \mathrm{nstar} }{2(b-a)}
% \left\|f'-\frac{f(b)-f(a)}{b-a}\right\|_1,\]
% 
% then the \(q\) output by this algorithm is guaranteed to satisfy
%
% \[\left\| \int_{a}^{b} f(x) dx - q \right\|_{1} \le \mathrm{abstol},\]
%
% provided the flag \(\mathrm{exceedbudget} = 0.\)
%
% And the upper bound of the cost is
%
% \[\sqrt{ \frac{\mathrm{nstar}* (b-a)^2 \mathrm{Var}(f')}{2 \times \mathrm{abstol}}}
% + 2 \times \mathrm{nstar} +4.\]
%
%
%% Examples
% *Example 1*

q = integral_g(@(x) x.^2)

% Integrate function x with default input parameter to make the error less
% than 1e-7.
%%
% *Example 2*

f = @(x) exp(-x.^2); q = integral_g(f,'a',1,'b',2,'nlo',100,'nhi',10000,...
    'abstol',1e-5,'nmax',1e7)

% Integrate function x^2 on [1,2] with lowest initial number of function 
% values 100 and highest initial number of function values 10000, absolute 
% error tolerance 1e-5 and cost budget 10000000.
%%
% *Example 3*

q = integral_g()

% Warning: Function f must be a function handle. Now GAIL is using 
% f(x)=exp(-100*(x-0.5)^2).
%% See Also
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/integral.html">integral</a>
% </html>
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/quad.html">quad</a>
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
% <a href="help_cubSobol_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_funmin_g.html">funmin_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell, Martha Razo, and Sunny Yun, "Reliable Adaptive
% Numerical Integration", 2015+, working.
%
% [2]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
%
% [3] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.2)
% [MATLAB Software], 2017. Available from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
