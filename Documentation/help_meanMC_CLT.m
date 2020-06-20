%% meanMC_CLT
% Monte Carlo method to estimate the mean of a random variable
%% Syntax
% sol = *MEANMC_CLT*(Y,absTol,relTol,alpha,nSig,inflate)
%
%% Description
%
% sol = *MEANMC_CLT*(Y,absTol,relTol,alpha,nSig,inflate) estimates the
%  mean, mu, of a random variable to within a specified error tolerance,
%  i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
%  1-alpha, where abstol is the absolute error tolerance.  The default
%  values are abstol=1e-2 and alpha=1%. Input Y is a function handle that
%  accepts a positive integer input n and returns an n x 1 vector of IID
%  instances of the random variable.
%
% This is a heuristic algorithm based on a Central Limit Theorem
%  approximation.
%   
% *Input Arguments*
%
% * Y --- the function or structure for generating n IID instances of a
%     random variable Y whose mean we want to estimate. Y is often defined
%     as a function of some random variable X with a simple distribution.
%     The input of Yrand should be the number of random variables n, the
%     output of Yrand should be n function values. For example, if Y = X.^2
%     where X is a standard uniform random variable, then one may define
%     Yrand = @(n) rand(n,1).^2.
%
% * absTol --- the absolute error tolerance, which should be
%  non-negative --- default = 1e-2
%
% * relTol --- the relative error tolerance, which should be
%  non-negative and no greater than 1 --- default = 0
%
% * alpha --- the uncertainty, which should be a small positive
%  percentage --- default = 1%
%
% * nSig --- the number of samples used to compute the sample variance
%  --- default = 1000
%
% * inflate --- the standard deviation inflation factor --- default = 1.2
%
% *Output Arguments*
%
% * Y --- the random generator
%
% * absTol --- the absolute error tolerance
%
% * relTol --- the relative error tolerance
%
% * alpha --- the uncertainty
%
% * mu --- the estimated mean of Y.
%
% * stddev --- sample standard deviation of the random variable
%
% * nSample --- total sample used.
%
% * time --- the time elapsed in seconds.
%
% * errBd --- the error bound.
%
%% Examples
%
%%
% *Example 1*
%
% Estimate the integral with integrand $f(x) = x_1 x_2$ in the interval
% $[0,1]^2$ with absolute tolerance $10^{-3}$ and relative tolerance $0$:

  [mu,out] = meanMC_CLT(@(n) prod(rand(n,2),2), 0.001);
  exact = 1/4;
  check = double(abs(exact - mu) < 2e-3)
  
%%
% *Example 2*
%
% Estimate the integral $f(x)=\exp(-x^2)$ in the interval $[0,1]$
% using $x$ as a control variate and relative error $10^{-3}$:

  f = @(x)[exp(-x.^2), x];
  YXn = @(n)f(rand(n,1));
  s = struct('Y',YXn,'nY',1,'trueMuCV',1/2);
  exact = erf(1)*sqrt(pi)/2;
  success = 0; runs = 1000; tol = 1e-3;
  for i=1:runs, success = success + double(abs(exact-meanMC_CLT(s,0,tol)) < tol*exact); end
  check = success >= 0.99 * runs  
  
%%
% *Example 3*
%
% Estimate the Keister's integration in dimension 1 with $a=1$, $\frac{1}{\sqrt 2}$ and
% using $\cos(x)$ as a control variate:

  normsqd = @(x) sum(x.*x,2);
  f = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)).* exp((1/2-a^2)*normt);
  f1 = @(x,a,d) f(normsqd(x),a,d);
  f2 = @(x)[f1(x,1,1),f1(x,1/sqrt(2),1),cos(x)];
  YXn = @(n)f2(randn(n,1));
  s = struct('Y',YXn,'nY',2,'trueMuCV',1/sqrt(exp(1)));
  [hmu,out] = meanMC_CLT(s,0,1e-3);
  exact = 1.380388447043143;
  check = double(abs(exact-hmu) < max(0,1e-3*abs(exact)))

%%
% *Example 4*
%
% Estimate the integral with integrand $f(x) = x_1^3 x_2^3 x_3^3$ in the
% interval $[0,1]^3$ with pure absolute error $10^{-3}$ using $x_1 x_2 x_3$ as
% a control variate:

  f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).*x(:,2).*x(:,3)];
  s = struct('Y',@(n)f(rand(n,3)),'nY',1,'trueMuCV',1/8);
  [hmu,out] = meanMC_CLT(s,1e-3,0);
  exact = 1/64;
  check = double(abs(exact-hmu) < max(1e-3,1e-3*abs(exact)))

%%
% *Example 5*
%
% Estimate the integrals with integrands $f_1(x) = x_1^3 x_2^3 x_3^3$ and 
% $f_2(x)= x_1^2 x_2^2 x_3^2- \frac{1}{27}+\frac{1}{64}$ in the interval $[0,1]^3$
% using $x_1 x_2 x_3$ and $x_1+x_2+x_3$ as control variates:

  f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, ...
            x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64, ...
            x(:,1).*x(:,2).*x(:,3), ...
            x(:,1)+x(:,2)+x(:,3)];
  s = struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5]);
  [hmu,out] = meanMC_CLT(s,1e-4,1e-3);
  exact = 1/64;
  check = double(abs(exact-hmu) < max(1e-4,1e-3*abs(exact)))
  
%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.