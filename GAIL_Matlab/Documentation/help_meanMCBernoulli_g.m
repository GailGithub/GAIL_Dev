%% meanMCBernoulli_g
% Monte Carlo method to estimate the mean of a Bernoulli random variable
% to within a specified error tolerance with guaranteed confidence level
% 1-alpha.
%% Syntax
%
% mu = *meanMCBernoulli_g*(Yrand)
%
% mu = *meanMCBernoulli_g*(Yrand,abstol,reltol,index,alpha,nmax)
%
% mu = *meanMCBernoulli_g*(Yrand,'abstol',abstol,'reltol',reltol,'errtype',errtype,
%   'alpha',alpha,'nmax',nmax)
%
% mu = *meanMCBernoulli_g*(Yrand,in_param)
%
% [mu, out_param] = *meanMCBernoulli_g*(Yrand,in_param)
%
%% Description
%
% mu = *meanMC_g*(Yrand) |estimates the mean of a Bernoulli random
% variable Y to within a specified error tolerance with guaranteed
% confidence level 99%. Input Yrand is a function handle that accepts a
% positive integer input n and returns an n x 1 vector of IID instances
% of the Bernoulli random variable Y.|
%
% mu = *meanMC_g*(Yrand,abstol,alpha,n_sigma,fudge,tbudget,nbudget,npcmax)
% | estimates the mean of a Bernoulli random variable Y to within an error
% tolerance with guaranteed confidence level 1-alpha using all ordered
% parsing inputs abstol, reltol, index, alpha and nmax.|
%
% mu =
% *meanMC_g*(Yrand,'abstol',abstol,'alpha',alpha,'fudge',fudge,'tbudget',
% tbudget,'nbudget',nbudget,'npcmax',npcmax,'checked',checked) |estimates
% the mean of a Bernoulli random variable Y to within a specified error
% tolerance with guaranteed confidence level 1-alpha. All the field-value
% pairs are optional and can be supplied in different order.|
%
% mu = *meanMC_g*(Yrand,in_param) |estimates the mean of a Bernoulli random
% variable Y to within a specified error tolerance in_param.abstol with
% guaranteed confidence level 1-alpha. If a field is not specified, the
% default value is used.|
%
% [mu, out_param] = *meanMC_g*(Yrand,in_param) |estimates the mean
% of a Bernoulli random variable Y to within a specified absolute error
% tolerance with the input parameters in_param and output parameters
% out_param.
%
% *Input Arguments*
%
% * Yrand --- |the function for generating IID instances of a Bernoulli
% random variable Y whose mean we want to estimate.|
%
% * p --- |the estimated mean of Y.|
%
% * in_param.abstol --- |the absolute error tolerance, default value is
% 1e-2.|
%
% * in_param.reltol --- |the relative error tolerance, default value is
% 1e-1.|
%
% * in_param.index --- |the error tolerance criterion, default value is
%   'abs'.|
%
% * in_param.alpha --- |the uncertainty, default value is 1%.|
%
% * in_param.nmax --- |the sample budget, default value is 1e8.|
%
% *Output Arguments*
%
% * out_param.n --- |the total sample used.|
%
% * out_param.time --- |the time elapsed.|
%
%% Examples
%
% Example 1:
% Calculate the mean of a bernoulli random variable with true p=0.55,with
% error tolerance 1e-3 and uncertainty 0.01.

in_param.abstol=1e-2; in_param.alpha = 0.01; p=1/90;Yrand=@(n) binornd(1,p,n,1);
p=meanMCBernoulli_g(Yrand,in_param)


%%
% Example 2:
% Using the same function as example 1, with the absolute error tolerance
% 1e-2.

p=meanMCBernoulli_g(Yrand,1e-2,1e-2,'abs')


%%
% Example 3:
% Using the sample function as example 1, with uncertainty 0.005.

p=meanMCBernoulli_g(Yrand,'index','rel','alpha',0.05)


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