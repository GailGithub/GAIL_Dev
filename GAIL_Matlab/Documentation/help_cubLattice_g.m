%% cubLattice_g
% Quasi-Monte Carlo method using rank-1 Lattices cubature
% over a d-dimensional region to integrate within a specified generalized error 
% tolerance with guarantees under Fourier coefficients cone decay assumptions.
%% Syntax
% [q,out_param] = *cubLattice_g*(f,d)
%
% q = *cubLattice_g*(f,d,abstol,reltol,measure,shift,mmin,mmax,fudge,transform,toltype,theta)
%
% q = *cubLattice_g*(f,d,'abstol',abstol,'reltol',reltol,'measure',measure,'shift',shift,'mmin',mmin,'mmax',mmax,'fudge',fudge,'transform',transform,'toltype',toltype,'theta',theta)
%
% q = *cubLattice_g*(f,d,in_param)
%% Description
%
% [q,out_param] = *cubLattice_g*(f,d) estimates the integral of f over the
%  d-dimensional region with an error guaranteed not to be greater than 
%  a specific generalized error tolerance, 
%  tolfun:=max(abstol,reltol*|integral(f)|). The generalized tolerance function can
%  aslo be cosen as tolfun:=theta*abstol+(1-theta)*reltol*|integral(f)| 
%  where theta is another input parameter. Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension of the hypercube,
%  and n is the number of points being evaluated simultaneously. The input d
%  is the dimension in which the function f is defined. Given the
%  construction of our Lattices, d must be a positive integer with 1<=d<=250.
% 
% q = *cubLattice_g*(f,d,abstol,reltol,measure,shift,mmin,mmax,fudge,transform,toltype,theta)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either.
% 
% q = *cubLattice_g*(f,d,'abstol',abstol,'reltol',reltol,'measure',measure,'shift',shift,'mmin',mmin,'mmax',mmax,'fudge',fudge,'transform',transform,'toltype',toltype,'theta',theta)
%  estimates the integral of f over a d-dimensional region. The answer
%  is given within the generalized error tolerance tolfun. All the field-value
%  pairs are optional and can be supplied with any order. If an input is not
%  specified, the default value is used.
% 
% q = *cubLattice_g*(f,d,in_param) estimates the integral of f over the
%  d-dimensional region. The answer is given within the generalized error 
%  tolerance tolfun.
% 
% *Input Arguments*
%
% * f --- the integrand whose input should be a matrix nxd where n is the
%  number of data points and d the dimension. By default it is the
%  quadratic function.
% 
% * d --- dimension of domain on which f is defined. d must be a positive
%  integer 1<=d<=250. By default it is 1.
% 
% * in_param.abstol --- the absolute error tolerance, abstol>0. By 
%  default it is 1e-4.
%
% * in_param.reltol --- the relative error tolerance, which should be
%  in [0,1]. Default value is 1e-1.
% 
% * in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%  measure of a uniformly distributed random variable in [0,1)^d
%  or normally distributed with covariance matrix I_d. By default it 
%  is 'uniform'. The only possible values are 'uniform' or 'normal'.
% 
% *Optional Input Arguments*
%
% * in_param.shift --- the Rank-1 lattices can be shifted to avoid the origin
%  or other particular points. By default we consider a uniformly [0,1)
%  random shift.
% 
% * in_param.mmin --- the minimum number of points to start is 2^mmin. The
%  cone condition on the Fourier coefficients decay requires a minimum
%  number of points to start. The advice is to consider at least mmin=10.
%  mmin needs to be a positive integer with mmin<=mmax. By default it is 10.
% 
% * in_param.mmax --- the maximum budget is 2^mmax. By construction of our
%  Lattices generator, mmax is a positive integer such that mmin<=mmax<=26.
%  The default value is 24.
% 
% * in_param.fudge --- the positive function multiplying the finite 
%  sum of Fast Fourier coefficients specified in the cone of functions.
%  For more information about this parameter, refer to the references.
%  By default it is @(x) 5*2.^-x.
% 
% * in_param.transform --- the algorithm is defined for continuous periodic functions. If the
%  input function f is not, there are 5 types of transform to periodize it
%  without modifying the result. By default it is Baker. The options:
%    'id' : no transformation. Choice by default.
%    'Baker' : Baker's transform or tent map in each coordinate. Preserving
%              only continuity but simple to compute.
%    'C0' : polynomial transformation only preserving continuity.
%    'C1' : polynomial transformation preserving the first derivative.
%    'C1sin' : Sidi transform with sinus preserving the first derivative.
%              This is in general a better option than 'C1'.
%
% * in_param.toltype --- this is the tolerance function. There are two
%  choices, 'max' (chosen by default) which takes
%  max(abstol,reltol*|integral(f)|) and 'comb' which is a linear combination
%  theta*abstol+(1-theta)*reltol*|integral(f)|. Theta is another 
%  parameter that can be specified (see below). For pure absolute error,
%  either choose 'max' and set reltol=0 or choose 'comb' and set
%  theta=1.
% 
% * in_param.theta --- this input is parametrizing the toltype 
%  'comb'. Thus, it is only afecting when the toltype
%  chosen is 'comb'. It stablishes the linear combination weight
%  between the absolute and relative tolerances
%  theta*abstol+(1-theta)*reltol*|integral(f)|. Note that for theta=1, 
%  we have pure absolute tolerance while for theta=0, we have pure 
%  relative tolerance. By default, theta=1.
%
% *Output Arguments*
%
% * q --- the estimated value of the integral.
% 
% * out_param.n --- number of points used when calling cubLattice_g for f.
% 
% * out_param.bound_err --- predicted bound on the error based on the cone
%  condition. If the function lies in the cone, the real error should be
%  smaller than this predicted error.
% 
% * out_param.time --- time elapsed in seconds when calling cubLattice_g for f.
%
% * out_param.exitflag --- this is a binary vector stating whether 
%  warning flags arise. These flags tell about which conditions make the
%  final result certainly not guaranteed. One flag is considered arisen
%  when its value is 1. The following list explains the flags in the
%  respective vector order:
%

%                    1    If reaching overbudget. It states whether
%                    the max budget is attained without reaching the
%                    guaranteed error tolerance.
%   
%                    2   If the function lies outside the cone. In
%                    this case, results are not guaranteed. Note that
%                    this parameter is computed on the transformed
%                    function, not the input function. For more
%                    information on the transforms, check the input
%                    parameter in_param.transfrom; for information about
%                    the cone definition, check the article mentioned
%                    below.
% 
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in
% dimension d with a prescribed generalized error tolerance. The Fourier
% coefficients of the integrand are assumed to be absolutely convergent. If
% the algorithm terminates without warning messages, the output is given
% with guarantees under the assumption that the integrand lies inside a
% cone of functions. The guarantee is based on the decay rate of the
% Fourier coefficients. For more details on how the cone is defined, please
% refer to the references below.
%
%% Examples
%
%%
% *Example 1*

% Estimate the integral with integrand f(x) = x1.*x2 in the interval
% [0,1)^2:

  f = @(x) prod(x,2); d = 2;
  q = cubLattice_g(f,d,1e-5,1e-1,'uniform','transform','C1sin')

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d = 3;
  q = cubLattice_g(f,d,1e-3,1e-3,'normal','transform','C1sin')

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); d = 2;
  q = cubLattice_g(f,d,1e-3,1e-1,'uniform','transform','C1')

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d = 1;
  q = cubLattice_g(f,d,1e-4,1e-1,'normal','fudge',@(m) 2.^-(2*m),...
      'transform','C1sin')

%%
% *Example 5*

% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the
% interval [0,1)^5 with pure absolute error 1e-5.

  f = @(x) 8*prod(x,2); d = 5;
  q = cubLattice_g(f,d,1e-5,0)
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
% Integration Based on Rank-1 Lattices (2014). Submitted for publication:
% arXiv:1411.1966.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
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
