%% meanMC_g
% |Monte Carlo method to estimate the mean of a random variable to within a
% specific absolute error tolerance with guaranteed uncertainty within alpha.|
%
%% Syntax
%
% mu = *meanMC_g*(Yrand)
%
% mu = *meanMC_g*(Yrand,abstol,alpha,n_sigma,fudge,tbudget,nbudget,npcmax,checked)
%
% mu = *meanMC_g*(Yrand,'abstol',abstol,'alpha',alpha,'fudge',fudge,'tbudget',
% tbudget,'nbudget',nbudget,'npcmax',npcmax,'checked',checked)
%
% mu = *meanMC_g*(Yrand,in_param)
%
% [mu, out_param] = *meanMC_g*(Yrand,in_param)
%
%% Description
%
% mu = *meanMC_g*(Yrand) |estimates the mean of a random variable Y to within
% a specified absolute error tolerance 1e-2 with guaranteed confidence
% level 99%. Input Yrand is a function handle that accepts a positive
% integer input n and returns an n x 1 vector of IID instances of the
% random variable Y.|
%
% mu = *meanMC_g*(Yrand,abstol,alpha,n_sigma,fudge,tbudget,nbudget,npcmax)
% | estimates the mean of a random variable Y to within an specified absolute
% error tolerance abstol with guaranteed confidence level 1-alpha. using
% all ordered parsing inputs abstol, n_sigma, fudge, tbudget, nbudget,
% npcmax and checked.|
%
% mu =
% *meanMC_g*(Yrand,'abstol',abstol,'alpha',alpha,'fudge',fudge,'tbudget',
% tbudget,'nbudget',nbudget,'npcmax',npcmax,'checked',checked) |estimates the mean of a random variable Y to within a
% specified absolute error tolerance abstol with guaranteed confidence
% level 1-alpha. All the field-value pairs are optional and can be
% supplied in different order.|
%
% mu = *meanMC_g*(Yrand,in_param) |estimates the mean of a random variable
% Y to within a specified absolute error tolerance in_param.abstol with
% guaranteed uncertainty within in_param.alpha. If a field is not
% specified, the default value is used.|
%
% [mu, out_param] = *meanMC_g*(Yrand,in_param) |estimates the mean of a
% random variable Y to within a specified absolute error tolerance with the
% given parameters in_param and output parameters out_param.|
%
% *Input Arguments*
%
% * Yrand --- |the function for generating IID instances of a random
%   variable Y whose mean we want to estimate. Y is often defined as a
%   function of some random variable X with a simple distribution.  For
%   example, if Y = X.^2 where X is a standard uniform random variable,
%   then one may define Yrand = @(n) rand(n,1).^2.|
%
% * mu --- |the estimated mean of Y.|
%
% * in_param.abstol --- |the absolute error tolerance, default value is
%   1e-2.|
%
% * in_param.alpha --- |the uncertainty, default value is 1%.|
%
% * in_param.n_sigma --- |initial sample size for estimating the sample
%   variance, the default value is 1e3.|
%
% * in_param.fudge --- |the standard deviation inflation factor, the
% default value is 1.1.|
%
% * in_param.tbudget --- |the time budget to do the two-stage
%   estimation, the default value is 100 seconds.|
%
% * in_param.nbudget --- |the sample budget to do the two-stage
%   estimation, the default value is 1e8.|
%
% * in_param.npcmax --- |number of elements in an array of optimal size to
%   calculate the mu, the default value is 1e6.|
%
% * in_param.checked --- |the value corresponding to parameter checking status.|
%
%                        0   not checked
%
%                        1   checked by cubMC_g
%
%                        2   checked by meanMC_g
%
% *Output Arguments*
%
% * out_param.time_n_sigma_predict --- |the estimated time to get n_sigma
%   samples of the random variable.|
%
% * out_param.n_left_predict --- |using the time left to predict the number
%   of samples left.|
%
% * out_param.nmax --- |the maximum sample budget to estimate mu, it comes
%   from both the sample budget and the time budget.|
%
% * out_param.var --- |the sample variance.|
%
% * out_param.exit --- |the state of program when exiting.|
%
%                        1  No enough samples to estimate the mean
%
%                        2  Initial try out time costs more than 10% of time budget
%
%                        3  The estimated time for estimating variance is bigger than 
%                           half of the time budget
%
% * out_param.kurtmax --- |the upper bound on modified kurtosis.|
%
% * out_param.n_mu --- |the sample needed to estimate the mu.|
%
% * out_param.n --- |the total sample needed to do the two stage estimation.|
%
% * out_param.time --- |the time elapsed.|
%
%% Guarantee
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
% Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G.
% W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to
% appear, arXiv:1208.4318 [math.ST]
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
% Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
% 1.3.0)" [MATLAB Software], 2014. Available from
% http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.