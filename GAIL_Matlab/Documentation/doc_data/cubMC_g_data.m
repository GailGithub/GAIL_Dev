%% Guarantee
%
% Error guarantee:
%
% Suppose the modified kurtosis, $\tilde{\kappa}$, of the integrand f
% satisfies the inequality:
%
% $$ \tilde{\kappa} \leq \frac{n_{\sigma}-3}{n_{\sigma}-1}+
% \left(\frac{\alpha n_\sigma}{1-\alpha}\right)\left(1-\frac{1}{C^2}\right)^2 =:
% \tilde{\kappa}_{\max} $$
%
% where $n_{\sigma}$ is the number of samples used to estimate the variance
% of f, C is the standard deviation inflation factor, and $\alpha$ is the
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
% integrands of variance no greater than $\sigma^2_{\max}$ and modified
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
%% Examples
%
% Example 1:
% Estimate the integral with integrand f(x) = sin(x) in the interval
% [1;2].

    f = @(x) sin(x);interval = [1;2]; Q = cubMC_g(f,interval,'uniform',1e-3)

%%
% Example 2:
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [0 0;1 1], where x is a vector x = [x1 x2].

    f = @(x) exp(-x(:,1).^2-x(:,2).^2);hyperbox = [0 0;1 1];
    Q = cubMC_g(f,hyperbox,'uniform',1e-3)
    
%%
% Example 3: 
% Estimate the integral with integrand $f(x) = 2^d\prod_{j=1}^d x_j+0.555$
% in the hyperbox [zeros(1,d);ones(1,d)], where x is a vector x = [x1 x2 ... xd].
%
    d=3;f=@(x) 2^d*prod(x,2)+0.555;hyperbox = [zeros(1,d);ones(1,d)];
    Q = cubMC_g(f,hyperbox,'uniform',1e-3)
