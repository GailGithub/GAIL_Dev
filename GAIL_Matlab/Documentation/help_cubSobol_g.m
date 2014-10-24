%% cubSobol_g
% |is a Quasi-Monte Carlo method using Sobol' cubature over the
% d-dimensional region to integrate within a specified absolute error 
% tolerance with guarantees under Walsh-Fourier coefficients cone decay assumptions.|
%% Syntax
% [q,out_param] = *cubSobol_g*(f,d)
%
% q = *cubSobol_g*(f,d,abstol,density,mmin,mmax,fudge)
%
% q = *cubSobol_g*(f,d,'abstol',abstol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge)
%
% q = *cubSobol_g*(f,d,in_param)
%% Description
%
% [q,out_param] = *cubSobol_g*(f,d) estimates the integral of f over the
%  d-dimensional region with an error guaranteed not to be greater than the
%  predefined error tolerance 1e-4. Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension of the hypercube,
%  and n is the number of points being evaluated simultaneously. The input d
%  is the dimension in which the function f is defined. Given the
%  construction of Sobol', d must be a positive integer with 1<=d<=1111.
%
% q = *cubSobol_g*(f,d,abstol,density,mmin,mmax,fudge)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the absolute error tolerance abstol. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail can not be specified either.
%
% q = *cubSobol_g*(f,d,'abstol',abstol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the absolute error tolerance abstol. All the field-value
%  pairs are optional and can be supplied with any order. If an input is not
%  specified, the default value is used.
%
% q = *cubSobol_g*(f,d,in_param) estimates the integral of f over the
%  d-dimensional region. The answer is given within the absolute error 
%  tolerance in_param.abstol.
% 
% *Input Arguments*
%
% * f --- |the integrand whose input should be a matrix nxd where n is the
%  number of data points and d the dimension. By default it is the
%  quadratic function.|
%
% * d --- |dimension of domain on which f is defined. d must be a positive
%  integer 1<=d<=1111. By default it is 1.|
%
% * in_param.abstol --- |the absolute error tolerance, abstol>0. By 
%  default it is 1e-4. |
%
% * in_param.density --- |for f(x), we can define x uniformly in [0,1)^d or
%  normally distributed with covariance matrix I_d. By default it is
%  'uniform'. The only possible values are 'uniform' or 'normal'.|
%
% * in_param.mmin --- |the minimum number of points to start is 2^mmin. The
%  cone condition on the Fourier coefficients decay requires a minimum
%  number of points to start. The advice is to consider at least mmin=10.
%  mmin needs to be a positive integer with mmin<=mmax. By default it is 10.|
%
% * in_param.mmax --- |the maximum budget is 2^mmax. By construction of the
%  Sobol' generator, mmax is a positive integer such that mmin<=mmax<=53.
%  The default value is 24.|
%
% * in_param.fudge --- |the positive function multiplying the finite 
%  sum of Fast Walsh coefficients specified in the cone of functions.
%  For more information about this parameter, refer to the references.
%  By default it is @(x) 5*2^-x.|
%
% *Output Arguments*
%
% * q --- |the estimated value of the integral.|
%
% * out_param.overbudget --- |boolean stating whether the max budget is
%  attained without reaching the guaranteed error tolerance. Output 1
%  means we have overrun our budget.|
%
% * out_param.n --- |number of points used when calling cubSobol_g for f.|
%
% * out_param.pred_err --- |predicted bound on the error based on the cone
%  condition. If the function lies in the cone, the real error should be
%  smaller than this predicted error.|
%
% * out_param.time --- |time elapsed in seconds when calling cubSobol_g for f.|
% 
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in [0,1)^d 
% with a prescribed absolute error tolerance. The Walsh-Fourier 
% coefficients of the integrand are assumed to be absolutely convergent.
% If the algorithm terminates without warning messages, the output is 
% given with guarantees under the assumption that the integrand lies inside
% a cone of functions. The guarantee is based on the decay rate of the 
% Walsh-Fourier coefficients. For more details on how the cone is defined, 
% please refer to the references below.
%
%% Examples
%
%%
% *Example 1*

% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:

  f=@(x) x(:,1).*x(:,2); d=2;
  q = cubSobol_g(f,d,1e-5,'uniform')

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f=@(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d=3;
  q = cubSobol_g(f,d,1e-3,'normal')

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:

  f=@(x) exp(-x(:,1).^2-x(:,2).^2); d=2;
  q = cubSobol_g(f,d,1e-3,'uniform')

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f=@(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d=1;
  q = cubSobol_g(f,d,1e-4,'normal','fudge',@(x) 2^-(2*x))
%% See Also
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama: Reliable adaptive
% cubature using digital sequences (2014). Submitted for publication.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% "GAIL: Guaranteed Automatic Integration Library (Version 2.0.0)"
% [MATLAB Software], 2014. Available from http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
