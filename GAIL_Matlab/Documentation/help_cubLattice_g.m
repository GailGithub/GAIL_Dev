%% cubLattice_g
% |is a Quasi-Monte Carlo method using rank-1 Lattices cubature
% over a d-dimensional region to integrate within a specified absolute error 
% tolerance with guarantees under Fourier coefficients cone decay assumptions.|
%% Syntax
% [q,out_param] = *cubLattice_g*(f,d)
%
% q = *cubLattice_g*(f,d,abstol,density,shift,mmin,mmax,fudge,transform)
%
% q = *cubLattice_g*(f,d,'abstol',abstol,'density',density,'shift',shift,'mmin',mmin,'mmax',mmax,'fudge',fudge,'transform',transform)
%
% q = *cubLattice_g*(f,d,in_param)
%% Description
%
% [q,out_param] = *cubLattice_g*(f,d) estimates the integral of f over the
%  d-dimensional region with an error guaranteed not to be greater than the
%  predefined error tolerance 1e-4. Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension of the hypercube,
%  and n is the number of points being evaluated simultaneously. The input d
%  is the dimension in which the function f is defined. Given the
%  construction of our Lattices, d must be a positive integer with 1<=d<=250.
% 
% q = *cubLattice_g*(f,d,abstol,density,shift,mmin,mmax,fudge,transform)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the absolute error tolerance abstol. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail can not be specified either.
% 
% q = *cubLattice_g*(f,d,'abstol',abstol,'density',density,'shift',shift,'mmin',mmin,'mmax',mmax,'fudge',fudge,'transform',transform)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the absolute error tolerance abstol. All the field-value
%  pairs are optional and can be supplied with any order. If an input is not
%  specified, the default value is used.
% 
% q = *cubLattice_g*(f,d,in_param) estimates the integral of f over the
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
%  integer 1<=d<=250. By default it is 1.|
% 
% * in_param.abstol --- |the absolute error tolerance, abstol>0. By 
%  default it is 1e-4.|
% 
% * in_param.density --- |for f(x)*mu(dx), we can define mu(dx) to be the
%  density function of a uniformly distributed random variable in [0,1)^d
%  or normally distributed with covariance matrix I_d. By default it 
%  is 'uniform'. The only possible values are 'uniform' or 'normal'.|
% 
% * in_param.shift --- |the Rank-1 lattices can be shifted to avoid the origin
%  or other particular points. By default we consider a uniformly [0,1)
%  random shift.|
% 
% * in_param.mmin --- |the minimum number of points to start is 2^mmin. The
%  cone condition on the Fourier coefficients decay requires a minimum
%  number of points to start. The advice is to consider at least mmin=10.
%  mmin needs to be a positive integer with mmin<=mmax. By default it is 10.|
% 
% * in_param.mmax --- |the maximum budget is 2^mmax. By construction of our
%  Lattices generator, mmax is a positive integer such that mmin<=mmax<=26.
%  The default value is 24.|
% 
% * in_param.fudge --- |the positive function multiplying the finite 
%  sum of Fast Fourier coefficients specified in the cone of functions.
%  For more information about this parameter, refer to the references.
%  By default it is @(x) 5*2^-x.|
% 
% * in_param.transform --- |the algorithm is defined for continuous periodic functions. If the
%  input function f is not, there are 5 types of transform to periodize it
%  without modifying the result. By default it is Baker. The options:
%    'id' : no transformation. Choice by default.
%    'Baker' : Baker's transform or tent map in each coordinate. Preserving
%              only continuity but simple to compute.
%    'C0' : polynomial transformation only preserving continuity.
%    'C1' : polynomial transformation preserving the first derivative.
%    'C1sin' : Sidi transform with sinus preserving the first derivative.
%              This is in general a better option than 'C1'.|
%
% *Output Arguments*
%
% * q --- |the estimated value of the integral.|
% 
% * out_param.overbudget --- |boolean stating whether the max budget is
%  attained without reaching the guaranteed error tolerance. Output 1
%  means we have overrun our budget.|
% 
% * out_param.n --- |number of points used when calling cubLattice_g for f.|
% 
% * out_param.pred_err --- |predicted bound on the error based on the cone
%  condition. If the function lies in the cone, the real error should be
%  smaller than this predicted error.|
% 
% * out_param.time --- |time elapsed in seconds when calling cubLattice_g for f.|
% 
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in [0,1)^d 
% with a prescribed absolute error tolerance. The Fourier coefficients of 
% the integrand are assumed to be absolutely convergent.
% If the algorithm terminates without warning messages, the output is 
% given with guarantees under the assumption that the integrand lies inside
% a cone of functions. The guarantee is based on the decay rate of the 
% Fourier coefficients. For more details on how the cone is defined, 
% please refer to the references below.
%
%% Examples
%
%%
% *Example 1*

% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:

  f=@(x) x(:,1).*x(:,2); d=2;
  q = cubLattice_g(f,d,1e-5,'uniform','transform','C1sin')

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f=@(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d=3;
  q = cubLattice_g(f,d,1e-3,'normal','transform','C1sin')

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:

  f=@(x) exp(-x(:,1).^2-x(:,2).^2); d=2;
  q = cubLattice_g(f,d,1e-3,'uniform','transform','C1')

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f=@(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d=1;
  q = cubLattice_g(f,d,1e-4,'normal','fudge',@(x) 2^-(2*x),'transform','C1sin')
%% See Also
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
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
% [1] Lluis Antoni Jimenez Rugama and Fred J. Hickernell: Adaptive Multidimensional
% Integration Based on Rank-1 Lattices (2014). Submitted for publication.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% "GAIL: Guaranteed Automatic Integration Library (Version 2.0.0)"
% [MATLAB Software], 2014. Available from http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
