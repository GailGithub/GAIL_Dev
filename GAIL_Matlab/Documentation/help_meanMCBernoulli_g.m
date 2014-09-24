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
% * in_param.errtype --- |the error type, default value is 'abs'.
%                       'abs'--- absolute error criterion
%                       'rel'--- relative error criterion
%                       'either'---absolute OR relative criterion|
% 
% * in_param.alpha --- |the uncertainty, default value is 1%.|
% 
% * in_param.nmax --- |the sample budget, default value is 1e9.|
% 
% *Output Arguments*
%
% * out_param.nabs --- |sample size needed to satisfy absolute error
%  tolerance|
%
% * out_param.nrel --- |sample size needed to satisfy relative error
%  tolerance |
%
% * out_param.n --- |the total sample used.|
%
% * out_param.tau --- |the iteration step.|
% 
% * out_param.time --- |the time elapsed.|
% 
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
% <a href="help_meanMC_g.html">meanMC_g</a>
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
%
