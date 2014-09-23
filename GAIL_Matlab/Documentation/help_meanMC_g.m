%% meanMC_g
% |MEANMC_G Monte Carlo method to estimate the mean of a random variable
%  to within a specified generalized error tolerance 
%  tolfun = max(abstol,reltol|mu|) with guaranteed confidence level 1-alpha.|
%% Syntax
% tmu = *meanMC_g*(Yrand)
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%
% tmu = *meanMC_g*(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha,
%
% [tmu, out_param] = *meanMC_g*(Yrand,in_param)
%% Description
%
% tmu = *meanMC_g*(Yrand) estimates the mean of a random variable Y to
%  within a specified generalized error tolerance tolfun =
%  max(abstol,reltol|mu|) with guaranteed confidence level 99%. Input
%  Yrand is a function handle that accepts a positive integer input n and
%  returns an n x 1 vector of IID instances of the random variable Y.
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%   estimates the mean of a random variable Y to within a specified
%   generalized error tolerance tolfun with guaranteed confidence
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
% * rand --- |the function for generating n IID instances of a random
%  ariable Y whose mean we want to estimate. Y is often defined as a
%  unction of some random variable X with a simple distribution. The
%  nput of Yrand should be the number of random variables n, the output
%  f Yrand should be n function values. For example, if Y = X.^2 where X
%  s a standard uniform random variable, then one may define Yrand =
%  (n) rand(n,1).^2.|
%
% * n_param.abstol --- |the absolute error tolerance, default value is 1e-2.|
%
% * n_param.reltol --- |the relative error tolerance, default value is 1e-1.|
%
% * n_param.alpha --- |the uncertainty, default value is 1%.|
%
% * n_param.fudge --- |standard deviation inflation factor, default value is
%  .2.|
%
% * n_param.nSig --- |initial sample size for estimating the sample
%  ariance, the default value is 1e3.|
%
% * n_param.n1 --- |initial sample size for estimating the sample
%  ean, the default value is 1e4.|
%
% * n_param.tbudget --- |the time budget to do the two-stage estimation,
%  he default value is 100 seconds.|
%
% * n_param.nbudget --- |the sample budget to do the two-stage estimation,
%  he default value is 1e9.|
%
% *Output Arguments*
%

%  mu--- the estimated mean of Y.|
%
% * ut_param.tau --- |the iteration step.|
%
% * ut_param.n --- |sample used in each iteration.|
%
% * ut_param.nmax --- |the maximum sample budget to estimate mu, it comes
%  rom both the sample budget and the time budget and sample has been
%  sed.|
%
% * ut_param.ntot --- |total sample used.|
%
% * ut_param.hmu --- |estimated mean in each iteration|
%
% * ut_param.tol --- |the tolerance for each iteration|
%
% * ut_param.var --- |the sample variance.|
%
% * ut_param.exit --- |the state of program when exiting.
%                  0   Success.
%                  1   Not enough samples to estimate the mean.
%                  2   Initial try out time costs more than 10% of time budget.
%                  3   The estimated time for estimating variance is bigger
%                      than half of the time budget.|
%
% * ut_param.kurtmax --- |the upper bound on modified kurtosis.|
%
% * ut_param.time --- |the time elapsed|
%
% * ut_param.checked --- |parameter checking status
%                     1  checked by meanMC_g|
%
%% Guarantee
%
% Error guarantee:
%
% Suppose the modified kurtosis, $\tilde{\kappa}$, of the random variable Y
% satisfies the inequality:
%
% $$\tilde{\kappa} \leq \frac{n_{\sigma}-3}{n_{\sigma}-1}+
% \left(\frac{\alpha n_\sigma}{1-\alpha}\right)\left(1-\frac{1}{C^2}\right)^2 =:
% \tilde{\kappa}_{\max}$$
%
% where $n_{\sigma}$ is the number of samples used to estimate the variance
% of Y, C is the standard deviation inflation factor, and $\alpha$ is the
% level of uncertainty. Then the answer $\hat{\mu}$ is guaranteed to
% satisfy the inequality:
%
% $$\mathrm{Pr}\left(|\mu-\hat{\mu}| \leq \varepsilon \right) \geq 1-\alpha$$
%
% where $\varepsilon$ is the absolute error tolerance.
%
% Cost upper bound guarantee:
%
% The probabilistic cost of the algorithm, with uncertainty $\beta$ , for
% random variables of variance no greater than $\sigma^2_{\max}$ and modified
% kurtosis no greater than $\tilde{\kappa}_{\max}$ is defined as
%
% $$N_{\mathrm{tot}}(\varepsilon,\alpha,\beta,\tilde{\kappa}_{\max},\sigma_{\max})
% := \sup_{\tilde{\kappa} \le \tilde{\kappa}_{\max}, \sigma \le
% \sigma_{\max} } \min\left\{N
% :\mathrm{Pr}[N_{\mathrm{tot}}(\varepsilon,\alpha,\tilde{\kappa}_{\max},\tilde{\kappa}_{\max}^{3/4})
% \le N] \ge 1-\beta  \right \}$$
%
% The total cost of this two stage algorithm has a probabilistic bound above
% by
%
% $$N_{\mathrm{tot}}(\varepsilon,\alpha, \beta, \tilde{\kappa}_{\max},
% \sigma_{\max}) \le N_{\mathrm{up}}(\varepsilon,\alpha, \beta,
% \tilde{\kappa}_{\max}, \sigma_{\max}) :=  n_{\sigma} +
% N_{\mu}(\varepsilon,\sigma_{\max}v(\tilde{\alpha},\beta,C),\tilde{\alpha},\tilde{\kappa}_{\max}^{3/4})
% $$
% with level of uncertainty $\beta$.
%
%% Examples
%
% Example 1:
% Calculate the mean of x^2 when x is uniformly distributed in [0 1], with
% the absolute error tolerance = 1e-2.

in_param.abstol=1e-2; in_param.alpha = 0.01; Yrand = @(n) rand(n,1).^2; 
mu = meanMC_g(Yrand,in_param)


%%
% Example 2:
% Using the same function as example 1, with the absolute error tolerance
% 1e-2.

mu = meanMC_g(Yrand,1e-2)


%%
% Example 3:
% Using the sample function as example 1, with the absolute error tolerance
% 1e-2 and uncertainty 0.01.

mu = meanMC_g(Yrand,'abstol',1e-3,'alpha',0.01)
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
%% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
% conservative fixed width confidence intervals via Monte Carlo sampling,
% Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
% Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
% arXiv:1208.4318 [math.ST]
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
% Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
% 1.3.0)" [MATLAB Software], 2014. Available from
% http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
