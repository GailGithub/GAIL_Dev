%% cubSobol_g
% |is a Quasi-Monte Carlo method using Sobol' cubature over the
% d-dimensional region to integrate within a specified generalized error
% tolerance with guarantees under Walsh-Fourier coefficients cone decay assumptions.|
%% Syntax
% [q,out_param] = *cubSobol_g*(f,d)
%
% q = *cubSobol_g*(f,d,abstol,reltol,density,mmin,mmax,fudge,errtype,theta)
%
% q = *cubSobol_g*(f,d,'abstol',abstol,'reltol',reltol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge,'errtype',errtype,'theta',theta)
%
% q = *cubSobol_g*(f,d,in_param)
%% Description
%
% [q,out_param] = *cubSobol_g*(f,d) estimates the integral of f over the
%  d-dimensional region with an error guaranteed not to be greater than 
%  a specific generalized error tolerance, 
%  tolfun:=max(abstol,reltol*|integral(f)|). The generalized tolerance function can
%  aslo be cosen as tolfun:=theta*abstol+(1-theta)*reltol*|integral(f)| 
%  where theta is another input parameter. Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension of the hypercube,
%  and n is the number of points being evaluated simultaneously. The input d
%  is the dimension in which the function f is defined. Given the
%  construction of Sobol', d must be a positive integer with 1<=d<=1111.
%
% q = *cubSobol_g*(f,d,abstol,reltol,density,mmin,mmax,fudge,errtype,theta)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either.
%
% q = *cubSobol_g*(f,d,'abstol',abstol,'reltol',reltol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge,'errtype',errtype,'theta',theta)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the generalized error tolerance tolfun. All the field-value
%  pairs are optional and can be supplied with any order. If an input is not
%  specified, the default value is used.
%
% q = *cubSobol_g*(f,d,in_param) estimates the integral of f over the
%  d-dimensional region. The answer is given within the generalized error 
%  tolerance tolfun.
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
%  default it is 1e-4.|
%
% * in_param.reltol --- |the relative error tolerance, which should be
%  in (0,1]. Default value is 1e-1.|
%
% * in_param.density --- |for f(x)*mu(dx), we can define mu(dx) to be the
%  density function of a uniformly distributed random variable in [0,1)^d
%  or normally distributed with covariance matrix I_d. By default it 
%  is 'uniform'. The only possible values are 'uniform' or 'normal'.|
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
%  By default it is @(x) 5*2.^-x.|
%
% * in_param.errtype --- |this is the tolerance function. There are two
%  choices, 'max' (chosen by default) which takes
%  max(abstol,reltol*|integral(f)|) and 'comb' which is a linear combination
%  theta*abstol+(1-theta)*reltol*|integral(f)|. Theta is another 
%  parameter that can be specified (see below).|
% 
% * in_param.theta --- |this input is parametrizing the errtype 
%  'comb'. Thus, it is only afecting when the errtype
%  chosen is 'comb'. It stablishes the linear combination weight
%  between the absolute and relative tolerances
%  theta*abstol+(1-theta)*reltol*|integral(f)|. Note that for theta=1, 
%  we have pure absolute tolerance while for theta=0, we have pure 
%  relative tolerance. By default, theta=1.|
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
% * out_param.outside_cone --- |boolean stating whether we did not meet
%  the necessary conditions for the integrand to be inside the cone. If
%  the value is true, the function is outside the cone. Otherwise, we do
%  not know. Note that this parameter is computed on the transformed
%  function, not the input function. For more information on the
%  transforms, check the input parameter in_param.transform.|
%
% * out_param.time --- |time elapsed in seconds when calling cubSobol_g for f.|
% 
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in dimension d 
% with a prescribed generalized error tolerance. The Walsh-Fourier 
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

  f = @(x) x(:,1).*x(:,2); d = 2;
  q = cubSobol_g(f,d,1e-5,1e-1,'uniform')

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d = 3;
  q = cubSobol_g(f,d,1e-3,1e-3,'normal')

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); d = 2;
  q = cubSobol_g(f,d,1e-3,1e-1,'uniform')

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d = 1;
  q = cubSobol_g(f,d,1e-4,1e-1,'normal','fudge',@(m) 2.^-(2*m))

%%
% *Example 5*

% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the interval
% [0,1)^5 with pure absolute error 1e-5.

  f = @(x) 8*prod(x,2); d = 5;
  q = cubSobol_g(f,d,1e-5,'errtype','comb','theta',1)
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
% cubature using digital sequences (2014). Submitted for publication:
% arXiv:1410.8615.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% "GAIL: Guaranteed Automatic Integration Library (Version 2.1)"
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
