%% cubSobol_g
% Quasi-Monte Carlo method using Sobol' cubature over the
% d-dimensional region to integrate within a specified generalized error
% tolerance with guarantees under Walsh-Fourier coefficients cone decay
% assumptions
%% Syntax
% [q,out_param] = *cubSobol_g*(f,hyperbox)
%
% q = *cubSobol_g*(f,hyperbox,measure,abstol,reltol)
%
% q = *cubSobol_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%
% q = *cubSobol_g*(f,hyperbox,in_param)
%% Description
%
% [q,out_param] = *cubSobol_g*(f,hyperbox) estimates the integral of f
%  over the d-dimensional region described by hyperbox, and with an error
%  guaranteed not to be greater than a specific generalized error tolerance,
%  tolfun:=max(abstol,reltol*| integral(f) |). Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension and n is the
%  number of points being evaluated simultaneously. 
%
% When measure is 'uniform', the input hyperbox is a 2 x d matrix, where
% the first row corresponds to the lower limits and the second row
% corresponds to the upper limits of the integral. When measure is 'uniform
% ball' or 'uniform sphere', the input hyperbox is a vector with d+1
% elements, where the first d values correspond to the center of the ball
% and the last value corresponds to the radius of the ball. For these last
% two measures, a user can optionally specify what transformation should be
% used in order to get a uniform distribution on a ball. When measure is
% 'uniform ball_box', the box-to-ball transformation, which gets a set of
% points uniformly distributed on a ball from a set of points uniformly
% distributed on a box, will be used. When measure is 'uniform
% ball_normal', the normal-to-ball transformation, which gets a set of
% points uniformly distributed on a ball from a set of points normally
% distributed on the space, will be used. Similarly, the measures 'uniform
% sphere_box' and 'uniform sphere_normal' can be used to specify the
% desired transformations. The default transformations are the box-to-ball
% and the box-to-sphere transformations, depending on the region of
% integration. Given the construction of Sobol' sequences, d must be a
% positive integer with 1 <= d<= 1111.
%
% q = *cubSobol_g*(f,hyperbox,measure,abstol,reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either. Inputs f and hyperbox
%  are required. The other optional inputs are in the correct order:
%  measure,abstol,reltol,mmin,mmax,and fudge.
%
% q = *cubSobol_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All the field-value
%  pairs are optional and can be supplied in any order. If an input is not
%  specified, the default value is used.
%
% q = *cubSobol_g*(f,hyperbox,in_param) estimates the integral of f over the
%  hyperbox. The answer is given within the generalized error tolerance tolfun.
%
% *Input Arguments*
%
% <html><ul type=square>
%  <li>f --- the integrand whose input should be a matrix n x d where n is
%  the number of data points and d the dimension, which cannot be
%  greater than 1111. By default f is f=@ x.^2.</li>
%  --- if using control variates, f needs to be a structure with two fields:
%  First field: 'func', need to be a function handle with n x (J+1)
%  dimension outputs, where J is the number of control variates.
%  First column is the output of target function, next J columns are
%  the outputs of control variates.
%  Second field: 'cv', need to be a 1 x J vector that stores the
%  exact means of control variates in the same order from
%  the function handle. For examples of how to use control variates,
%  please check Example 7 below.</li>
%  </ul>
% </html>
%
% * hyperbox --- the integration region defined by its bounds. When measure
%  is 'uniform' or 'normal', hyperbox must be a 2 x d matrix, where the
%  first row corresponds to the lower limits and the second row corresponds
%  to the upper limits of the integral. When measure is 'uniform ball'
%  or 'uniform sphere', the input hyperbox is a vector with d+1 elements,
%  where the first d values correspond to the center of the ball and the
%  last value corresponds to the radius of the ball. The default value
%  is [0;1].
%
% * in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%  measure of a uniformly distributed random variable in the hyperbox
%  or normally distributed with covariance matrix I_d. The possible
%  values are 'uniform', 'normal', 'uniform ball', 'uniform ball_box',
%  'uniform ball_normal', 'uniform sphere', 'uniform sphere_box' and
%  'uniform sphere_normal'. For 'uniform', the hyperbox must be a
%  finite volume, for 'normal', the hyperbox can only be defined as
%  (-Inf,Inf)^d and, for 'uniform ball' or 'uniform sphere', hyperbox
%  must have finite values for the coordinates of the center and a
%  finite positive value for the radius. By default it is 'uniform'.
%
% * in_param.abstol --- the absolute error tolerance, abstol>=0. By
%  default it is 1e-4. For pure absolute tolerance, set in_param.reltol
%  = 0.
%
% * in_param.reltol --- the relative error tolerance, which should be
%  in [0,1]. Default value is 1e-2. For pure absolute tolerance, set
%  in_param.abstol = 0.
%
% *Optional Input Arguments*
%
% * in_param.mmin --- the minimum number of points to start is 2^mmin.
%  The cone condition on the Fourier coefficients decay requires a
%  minimum number of points to start. The advice is to consider at least
%  mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%  default it is 10.
%
% * in_param.mmax --- the maximum budget is 2^mmax. By construction of
%  the Sobol' generator, mmax is a positive integer such that
%  mmin<=mmax<=53. The default value is 24.
%
% * in_param.fudge --- the positive function multiplying the finite
%  sum of Fast Walsh Fourier coefficients specified in the cone of functions.
%  This input is a function handle. The fudge should accept an array of
%  nonnegative integers being evaluated simultaneously. For more
%  technical information about this parameter, refer to the references.
%  By default it is @(m) 5*2.^-m.
%
% *Output Arguments*
%
% * q --- the estimated value of the integral.
%
% * out_param.d --- dimension over which the algorithm integrated.
%
% * out_param.n --- number of Sobol' points used for computing the
%  integral of f.
%
% * out_param.bound_err --- predicted bound on the error based on the cone
%  condition. If the function lies in the cone, the real error will be
%  smaller than generalized tolerance.
%
% * out_param.time --- time elapsed in seconds when calling cubSobol_g.
%
% * out_param.beta --- the value of beta when using control variates
%     as in f-(h-Ih)beta, if using 'betaUpdate' option, beta is a vector
%     storing value of each iteration.
%
% * y --- fast transform coefficients of the input function.
%
% * kappanumap --- wavenumber mapping used in the error bound.
%
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- this is a binary vector stating whether
%  warning flags arise. These flags tell about which conditions make the
%  final result certainly not guaranteed. One flag is considered arisen
%  when its value is 1. The following list explains the flags in the
%  respective vector order:</li>
%   <ul type="circle">
%                <li>1 : If reaching overbudget. It states whether
%                the max budget is attained without reaching the
%                guaranteed error tolerance.</li>
%                <li>2 : If the function lies outside the cone. In
%                this case, results are not guaranteed. For more
%                information about the cone definition, check the
%                article mentioned below.</li>
%   </ul>
%  </ul>
% </html>
%
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in $[0,1]^d$
% with a prescribed generalized error tolerance. The Walsh-Fourier
% coefficients of the integrand are assumed to be absolutely convergent. If
% the algorithm terminates without warning messages, the output is given
% with guarantees under the assumption that the integrand lies inside a
% cone of functions. The guarantee is based on the decay rate of the
% Walsh-Fourier coefficients. For integration over domains other than
% $[0,1]^d$, this cone condition applies to $f \circ \psi$ (the
% composition of the functions) where $\psi$ is the transformation
% function for $[0,1]^d$ to the desired region. For more details on how the
% cone is defined, please refer to the references below.
%
%% Examples
%
%%
% *Example 1*
%
% Estimate the integral with integrand $f(x) = x_1 x_2$ in the hyperbox $[0,1]^2$:

  f = @(x) prod(x,2); hyperbox = [zeros(1,2); ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
  check = double(abs(exactsol-q) < 1e-5)

%%
% *Example 2*
%
% Estimate the integral with integrand $f(x) = x_1^2  x_2^2 x_3^2$
% in the hyperbox $R^3$ where $x_1$, $x_2$ and $x_3$ are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; 
  hyperbox = [-inf(1,3);inf(1,3)]; exactsol = 1;
  q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3, 'fudge',@(m) 10*2.^-m); 
  check = double(abs(exactsol-q) < max(1e-3,1e-3*abs(exactsol)))

%%
% *Example 3*
%
% Estimate the integral with integrand $f(x) = exp(-x_1^2-x_2^2)$ in the
% hyperbox $[-1,2]^2$:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2); 2*ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-3,1e-2); 
  exactsol = 1/9*(sqrt(pi)/2*(erf(2)+erf(1)))^2;
  check = double(abs(exactsol-q) < max(1e-3,1e-2*abs(exactsol)))

%%
% *Example 4*
%
% Estimate the price of an European call with $S_0=100$, $K=100$, $r=\sigma^2/2$,
% $\sigma=0.05$, and $T=1$.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); 
  hyperbox = [-inf(1,1);inf(1,1)];
  q = cubSobol_g(f,hyperbox,'normal',1e-4,1e-2); 
  price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
  check = double(abs(price-q) < max(1e-4,1e-2*abs(price)))

%%
% *Example 5*
%
% Estimate the integral with integrand $f(x) = 8 x_1 x_2 x_3 x_4 x_5$ in the interval
% $[0,1)^5$ with pure absolute error $10^{-5}$.

  f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
  check = double(abs(exactsol-q) < 1e-5)

%%
% *Example 6*
%
% Estimate the integral with integrand $f(x) = x_1^2+x_2^2$ over the disk
% with center $(0,0)$ and radius $1$ with pure absolute error
% $10^{-5}$, where $x = [x_1, x_2]$ is a vector.

  f = @(x) x(:,1).^2+x(:,2).^2; hyperbox = [0,0,1]; exactsol = pi/2;
  q = cubSobol_g(f,hyperbox,'uniform ball','abstol',1e-4,'reltol',0); 
  check = double(abs(exactsol-q) < 1e-4)

%%
% *Example 7*
%
% Estimate the integral with integrand $f(x) = 10 x_1 - 5 x_2^2 + x_3^3$
% in the interval $[0,2)^3$ with pure absolute error $10^{-5}$ using
% two control variates $h_1(x) = x_1$ and $h_2(x) = x_2^2$.

  g.func = @(x) [10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
  g.cv = [1,4/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
  q = cubSobol_g(g,hyperbox,'uniform',1e-6,0); exactsol = 16/3;
  check = double(abs(exactsol-q) < 1e-6)


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
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <p><a href="cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama "Reliable
% adaptive cubature using digital sequences", Monte Carlo and Quasi-Monte
% Carlo Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools and D.
% Nuyens, eds.), Springer Proceedings in Mathematics and Statistics, vol.
% 163, Springer-Verlag, Berlin, 2016, arXiv:1410.8615 [math.NA], pp.
% 367-383.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice and Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% [6] Fang, K.-T., & Wang, Y. (1994). Number-theoretic Methods in
% Statistics. London, UK: CHAPMAN & HALL
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
