%% meanMCBernoulli_g
% |Monte Carlo method to estimate the mean of a Bernoulli
% random variable to within a specified error tolerance with guaranteed
% confidence level 1-alpha.|
%% Syntax
% pHat = *meanMCBernoulli_g*(Yrand)
%
% pHat = *meanMCBernoulli_g*(Yrand,abstol,reltol,errtype,alpha,nmax)
%
% pHat = *meanMCBernoulli_g*(Yrand,'abstol',abstol,'reltol',reltol,
%
% [pHat, out_param] = *meanMCBernoulli_g*(Yrand,in_param)
%% Description
%
% pHat = *meanMCBernoulli_g*(Yrand) estimates the mean of a Bernoulli random
%  variable Y to within a specified error tolerance with guaranteed
%  confidence level 99%. Input Yrand is a function handle that accepts a
%  positive integer input n and returns a n x 1 vector of IID instances
%  of the Bernoulli random variable Y.
% 
% pHat = *meanMCBernoulli_g*(Yrand,abstol,reltol,errtype,alpha,nmax) estimates
%  the mean of a Bernoulli random variable Y to within a specified error
%  tolerance with guaranteed confidence level 1-alpha using all ordered
%  parsing inputs abstol, reltol, errtype, alpha and nmax.
% 
% pHat = *meanMCBernoulli_g*(Yrand,'abstol',abstol,'reltol',reltol,
%  'errtype',errtype,'alpha',alpha,'nmax',nmax) estimates the mean of a
%  Bernoulli random variable Y to within a specified error tolerance with
%  guaranteed confidence level 1-alpha. All the field-value pairs are
%  optional and can be supplied in different order.
% 
% [pHat, out_param] = *meanMCBernoulli_g*(Yrand,in_param) estimates the mean
%  of a Bernoulli random variable Y to within a specified error tolerance
%  with the given parameters in_param and produce the estimated mean pHat
%  and output parameters out_param.
% 
% *Input Arguments*
%
% * Yrand --- |the function for generating IID instances of a Bernoulli
%            random variable Y whose mean we want to estimate.|
%
% * pHat --- |the estimated mean of Y.|
%
% * in_param.abstol --- |the absolute error tolerance, default value is 1e-2.|
%
% * in_param.reltol --- |the relative error tolerance, default value is 1e-1.|
%
% * in_param.errtype --- |the error type, default value is 'abs'.|
%   
%                       'abs'--- absolute error criterion
%   
%                       'rel'--- relative error criterion
%   
%                       'either'---absolute OR relative criterion
% 
% * in_param.alpha --- |the uncertainty, default value is 1%.|
% 
% * in_param.nmax --- |the sample budget, default value is 1e9.|
% 
% *Output Arguments*
%
% * out_param.nabs --- |sample size needed to satisfy absolute error 
%                     tolerance.|
%
% * out_param.nrel --- |sample size needed to satisfy relative error
%                     tolerance.|
%
% * out_param.n --- |the total sample used.|
%
% * out_param.tau --- |the iteration step.|
% 
% * out_param.time --- |the time elapsed.|
% 
%  Guarantee
%
% Case 1: errtype = 'abs'
%
% If the sample size is calculated according Hoeffding's inequality, which
% equals to ceil(log(2/out_param.alpha)/(2*out_param.abstol^2)), then the
% following inequality must be satisfied:
%
% Pr(|p-pHat| <= abstol) >= 1-alpha.
% 
% Here p is the true mean of Yrand, and pHat is the output of
% MEANMCBERNOULLI_G with errtype = 'abs'
%
% Also, the cost is deterministic and bounded.
% 
% Case 2: errtype = 'rel'
%
% If the algorithm terminated without any warning messages, the estimated
% mean pHat would satisfy the following inequality:
%
% Pr(|p-pHat| <= abstol*p) >= 1-alpha.
%
% Here p is the true mean of Y, and pHat is the output of MEANMCBERNOULLI_G
% with errtype = 'rel'.
% 
% Additionally, the cost of the algorithm would be bounded by N_up, which is
% defined in terms of the true mean p, uncertainty alpha and relative
% tolerance reltol. For details, please refer to the paper.
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
