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
