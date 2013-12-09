%% cubMC_g
% |Monte Carlo method to evaluate a multidimentional integral to
% within a specified absolute error tolerance with guaranteed uncertainy
% within alpha.|
%
% |The guarantee holds if the modified kurtosis is less than the kmax,
% i.e.:|
%
% $$\tilde{k} \leq \frac{n_{\sigma}-3}{n_{\sigma}-1}+
% (\frac{\alpha n_\sigma}{1-\alpha})(1-\frac{1}{c^2})^2: = k_{max}$$
%
%% Syntax
% [Q,out_param] = *cubMC_g*(f)
%
% Q = *cubMC_g*(f,interval,measure,abstol,alpha,n_sigma,fudge)
%
% Q = *cubMC_g*(f,interval,'measure','uniform','abstol',abstol,'alpha',alpha,
%              'n_sigma',n_sigma,fudge',fudge)
%
% Q = *cubMC_g*(f,interval,,in_param)
%
%% Description
% [Q,out_param] = *cubMC_g*(f,interval) |estimates the integral with
% integrand f to within the absolute error tolerance 1e-2 and with
% guaranteed uncertainty alpha within 1%. Input f a function handle. The
% function Y=f(X) should accept a vector argument X and return a vector
% result Y, the integrand evaluated at each element of X.|
%
% Q = *cubMC_g*(f,interval,measure,abstol,alpha,n_sigma,fudge) |estimates the
% integral with integrand f to within an absolute error tolerance abstol
% with guaranteed uncertainty within alpha using ordered parameter input
% interval, measure, tolerence, uncertainty, n_sigma and fudge factor.|
%
% Q = *cubMC_g*(f,interval,'measure','uniform','abstol',abstol,'alpha',alpha,
% 'n_sigma',n_sigma,fudge',fudge) |estimates the integral with integrand f
% to within an absolute error tolerance abstol with guaranteed uncertainty
% within alpha. All the field-value pairs are optional and can be supplied
% in different order.|
%
% Q = *cubMC_g*(f,interval,in_param) |estimates the integral with integrand f
% to within an absolute error tolerance in_param.abstol with guaranteed
% uncertainty within in_param.alpha. If a field is not specified, the
% default value is used.|
%
% *Input Arguments*
%
% * f --- |the integrand.|
%
% * interval --- |the integration interval.|
%
% * in_param.measure --- |the measure for generating the random variable,
%   the default is uniform.|
%
% * in_param.abstol --- |the absolute error tolerance, default value is 1e-2.|
%
% * in_param.alpha --- |the uncertainty, default value is 1%.|
%
% * in_param.n_sigma --- |initial sample size for estimating the sample
%                         variance, the default value is 1e3.|
%
% * in_param.fudge --- |the standard deviation inflation factor, the
%                       default value is 1.1.|
%
% *Output Arguments*
%
% * Q --- |the estimated value of the the integration.|
%
% * out_param_time_n_sigma_predict --- |the estimated time to get n_sigma
%                                       samples of the random variable.|
%
% * out_param.n_left_predict --- |using the time left to predict the number
%                                 of samples left.|
%
% * out_param.nmax --- |the maximum sample budget to estimate mu, it comes
%                       from both the sample budget and the time budget.|
%
% * out_param.var --- |the sample variance.|
%
% * out_param.kurtmax --- |the upper bound on modified kurtosis.|
%
% * out_param.time --- |the time eclipsed.|
%
% * out_param.n_mu --- |the sample size that needed to estimate the mu.|
%
% * out_param.n --- |the total sample size needed to do the two stage
% algorithm.|
%
% * out_param.exit --- |the state of program when exiting.|
%
%                         0   success
%
%                         1   No enough samples to estimate the mean
%
%                         2   Initial try out time costs more than
%                             10% of time budget
%
%                         3   The estimated time for estimating variance 
%                             is bigger than half of the time budget
%
%                         10  Interval does not contain numbers
%
%                         11  Interval not 2 x d
%
%                         12  Interval is only a point in one direction
%
%                         13  Interval is infinite when measure is uniform
%
%                         14  Interval is not doubly infinite when measure
%                             is normal
%
%% Examples
% Example 1:
% Estimate the integral with integrand $f(x) = x.^2$ in the interval [0,1].

    f = @(x) x.^2;interval = [0;1]; Q = cubMC_g(f,interval,'abstol',1e-2)
%%
% Example 2:
% Estimate the integral with integrand $f(x) = \exp(x)$ in the interval
% [1,2].

    f = @(x) exp(x);interval = [1;2]; Q = cubMC_g(f,interval)

%%
% Example 3:
% Estimate the integral with integrand $f(x) = \sin(x)$ in the interval
% [1,2].

    f = @(x) sin(x);interval = [1;2]; Q = cubMC_g(f,interval,'uniform',1e-3)

%%
% Example 4:
% Estimate the integral with integrand $f(x) = \exp(-x1^2-x2^2)$ in the
% interval [0 0;1 1],where x is a vector x = [x1 x2].

    f = @(x) exp(-x(:,1).^2-x(:,2).^2);interval = [0 0;1 1];
    Q = cubMC_g(f,interval,'uniform',1e-3)
    
%%
% Example 5: 
% Estimate the integral with integrand $f(x) = 2^d\prod_{j=1}^d x_j+0.555$
% in the interval [zeros(1,d);ones(1,d)], where x is a vector x = [x1 x2 ... xd].
%
    d=3;f=@(x) 2^d*prod(x,2)+0.555;interval = [zeros(1,d);ones(1,d)];
    Q = cubMC_g(f,interval,'uniform',1e-3)
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
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
%% Reference
%   [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G.
%   W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to
%   appear, arXiv:1208.4318 [math.ST]