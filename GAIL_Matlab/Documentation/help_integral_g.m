%% integral_g
% 1-D guaranteed function integration using trapezoidal rule
%
%% Syntax
% q = *integral_g*(f)
%
% q = *integral_g*(f,abstol,ninit,nmax)
%   
% q = *integral_g*(f,'abstol',abstol,'ninit',ninit,'nmax',nmax)
%
% q = *integral_g*(f,in_param)
%
% [q, out_param] = *integral_g*(f,...)
%
%% Description
%
% q = *integral_g*(f) |computes q, the definite integral of function f
% on the interval [0,1] by trapezoidal rule with in a guaranteed absolute
% error of 1e-6. Default starting number of sample points taken is 52 and 
% default cost budget is 1e7. Input f is a function handle. The function 
% y = f(x) should accept a vector argument x and return a vector result 
% y, the integrand evaluated at each element of x.|
%
% q = *integral_g*(f,in_param) |computes q, the definite integral of 
% function f by trapezoidal rule within a guaranteed absolute error
% in_param.abstol, starting number of points in_param.ninit, and cost 
% budget in_param.nmax. If a field is not specified, the default value
% is used.|
% 
% q = *integral_g*(f,'abstol',abstol,'ninit',ninit,'nmax',nmax) |computes 
% q, the definite integral of function f by trapezoidal rule within a 
% guaranteed absolute error tolerance abstol, starting number of points 
% ninit, and cost budget nmax. All three field-value pairs are optional 
% and can be supplied.|
%
% q = *integral_g*(f,abstol,ninit, nmax) |computes q, the defintie 
% integral of function f by trapezoidal rule with the ordered input 
% parameters, guaranteed absolute error tolerance abstol, starting number
% of points ninit, and cost budget nmax.|
%
% [q, out_param] = *integral_g*(f,...) |returns the approximated 
% integration q and output structure out_param, which includes the
% fields in_param.|
%
% *Input Arguments*
%   
% * in_param.abstol --- |absolute error tolerance required by user.|
% 
% * in_param.ninit --- |starting number of points specified by user.|
%
% * in_param.nmax --- |cost budget required by user.|
%
% *Output Arguments*
%
% * out_param.exceedbudget --- |it is true if the algorithm tries to use 
% more points than cost budget, false otherwise.|
%
% * out_param.ninit --- |the minimum requirement of number of points to 
% start|
% 
% * out_param.tauchange --- |it is true if the cone constant has been
% changed, false otherwise. See [1] for details.  If true, you may wish 
% to change the input in_param.ninit to a larger number.|
% 
% * out_param.npoints --- |number of points we need to 
% reach the guaranteed absolute error tolerance abstol.|
%
% * out_param.errest --- |approximation error defined as the differences
% between the true value and the approximated value of the integral.|
% 
%% Examples
%  Example 1

f = @(x) x; [q, out_param] = integral_g(f)

% Integrate function x with default input parameter to make the error
% less than 1e-7.

%%
% Example 2

[q, out_param] = integral_g(@(x) 3*x.^2,'ninit',52,'nmax',1e7,'abstol',1e-8)

% Integrate function x^2 with starting number of points 52, cost budget 10000000
% and error toerence 1e-8


%% See Also
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
%% Reference
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
%      The complexity of guaranteed automatic algorithms: Cones, not
%      balls, Journal of Complexity 2013, to appear, DOI
%      10.1016/j.jco.2013.09.002.