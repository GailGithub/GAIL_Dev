%% meanMC_g
% |MEANMC_G Monte Carlo method to estimate the mean of a random variable.|
%% Syntax
% tmu = *meanMC_g*(Yrand)
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%
% tmu = *meanMC_g*(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha,'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',tbudget,'nbudget',nbudget)
%
% [tmu, out_param] = *meanMC_g*(Yrand,in_param)
%% Description
%
% tmu = *meanMC_g*(Yrand) estimates the mean, mu, of a random variable Y to
%  within a specified generalized error tolerance, 
%  tolfun:=max(abstol,reltol*|mu|), i.e., |mu - tmu| <= tolfun with
%  probability at least 1-alpha, where abstol is the absolute error
%  tolerance, and reltol is the relative error tolerance. Usually the
%  reltol determines the accuracy of the estimation, however, if the |mu|
%  is rather small, the abstol determines the accuracy of the estimation.
%  The default values are abstol=1e-2, reltol=1e-1, and alpha=1%. Input
%  Yrand is a function handle that accepts a positive integer input n and
%  returns an n x 1 vector of IID instances of the random variable Y.
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%  estimates the mean of a random variable Y to within a specified
%  generalized error tolerance tolfun with guaranteed confidence
%  level 1-alpha using all ordered parsing inputs abstol, reltol, alpha,
%  fudge, nSig, n1, tbudget, nbudget.
%
% tmu = *meanMC_g*(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha,
%  'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',tbudget,'nbudget', nbudget)
%  estimates the mean of a random variable Y to within a specified
%  generalized error tolerance tolfun with guaranteed confidence level
%  1-alpha. All the field-value pairs are optional and can be supplied in
%  different order, if a field is not supplied, the default value is used.
%
% [tmu, out_param] = *meanMC_g*(Yrand,in_param) estimates the mean of a
%  random variable Y to within a specified generalized error tolerance
%  tolfun with the given parameters in_param and produce the estimated
%  mean tmu and output parameters out_param. If a field is not specified,
%  the default value is used.
%
% *Input Arguments*
%
% * Yrand --- |the function for generating n IID instances of a random
%  variable Y whose mean we want to estimate. Y is often defined as a
%  function of some random variable X with a simple distribution. The
%  input of Yrand should be the number of random variables n, the output
%  of Yrand should be n function values. For example, if Y = X.^2 where X
%  is a standard uniform random variable, then one may define Yrand =
%  @(n) rand(n,1).^2.|
%
% * in_param.abstol --- |the absolute error tolerance, which should be
%  positive, default value is 1e-2.|
%
% * in_param.reltol --- |the relative error tolerance, which should be
%  between 0 and 1, default value is 1e-1.|
%
% * in_param.alpha --- |the uncertainty, which should be a small positive
%  percentage. default value is 1%.|
%
% * in_param.fudge --- |standard deviation inflation factor, which should
%  be larger than 1, default value is 1.2.|
%
% * in_param.nSig --- |initial sample size for estimating the sample
%  variance, which should be a moderate large integer at least 30, the
%  default value is 1e4.|
%
% * in_param.n1 --- |initial sample size for estimating the sample mean,
%  which should be a moderate large positive integer at least 30, the
%  default value is 1e4.|
%
% * in_param.tbudget --- |the time budget in seconds to do the two-stage
%  estimation, which should be positive, the default value is 100 seconds.|
%
% * in_param.nbudget --- |the sample budget to do the two-stage
%  estimation, which should be a large positive integer, the default
%  value is 1e9.|
%
% *Output Arguments*
%
% * tmu --- |the estimated mean of Y.|
%
% * out_param.tau --- |the iteration step.|
%
% * out_param.n --- |the sample size used in each iteration.|
%
% * out_param.nmax --- |the maximum sample budget to estimate mu. It was
%   calculated by the sample left and time left.|
%
% * out_param.ntot --- |total sample used.|
%
% * out_param.hmu --- |estimated mean in each iteration.|
%
% * out_param.tol --- |the reliable upper bound on error for each iteration.|
%
% * out_param.var --- |the sample variance.|
%
% * out_param.exit --- |the state of program when exiting.|
%    
%                   0   Success
%   
%                   1   Not enough samples to estimate the mean
%
% * out_param.kurtmax --- |the upper bound on modified kurtosis.|
%
% * out_param.time --- |the time elapsed in seconds.|
%
% * out_param.flag --- |parameter checking status
%   
%                        1  checked by meanMC_g|
%
%%  Guarantee
% This algorithm attempts to calculate the mean of a random variable to a
% prescribed error tolerance with guaranteed confidence level 1-alpha. If
% the algorithm terminated without showing any warning messages and provide
% an answer tmu, then the follow inequality would be satisfied:
% 
% Pr(|mu-tmu| <= max(abstol,reltol*|mu|)) >= 1-alpha
%
% where abstol is the absolute error tolerance and reltol is the relative
% error tolerance, if the true mean mu is rather small as well as the
% reltol, then the abstol would be satisfied, and vice versa. 
%
% The cost of the algorithm is also bounded above by N_up, which is defined
% in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And the
% following inequality holds:
% 
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
%
%% Examples
%
%%
% *Example 1*

% Calculate the mean of x^2 when x is uniformly distributed in
% [0 1], with the relative error tolerance = 1e-3 and uncertainty 5%.

  in_param.reltol=0; in_param.abstol = 1e-3;
  in_param.alpha = 0.05; Yrand=@(n) rand(n,1).^2;
  tmu = meanMC_g(Yrand,in_param)

%%
% *Example 2*

% Calculate the mean of exp(x) when x is uniformly distributed in
% [0 1], with the absolute error tolerance 1e-3.

  tmu = meanMC_g(@(n)exp(rand(n,1)),1e-3,0)

%%
% *Example 3*

% Calculate the mean of sin(x) when x is uniformly distributed in
% [0 1], with the relative error tolerance 1e-2 and uncertainty 0.05.

  tmu = meanMC_g(@(n)cos(rand(n,1)),'reltol',1e-2,'abstol',0,'alpha',0.05)
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
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMCBernoulli_g.html">meanMCBernoulli_g</a>
% </html>
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
%% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
% conservative fixed width confidence intervals via Monte Carlo sampling,
% Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
% Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
% arXiv:1208.4318 [math.ST]
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou, "GAIL:
% Guaranteed Automatic Integration Library (Version 2.0)" [MATLAB
% Software], 2014. Available from http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
