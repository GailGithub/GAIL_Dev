%% meanMCBer_g
% |Monte Carlo method to estimate the mean of a Bernoulli random
% variable to within a specified absolute error tolerance with guaranteed
% confidence level 1-alpha.|
%% Syntax
% pHat = *meanMCBer_g*(Yrand)
%
% pHat = *meanMCBer_g*(Yrand,abstol,alpha,nmax)
%
% pHat = *meanMCBer_g*(Yrand,'abstol',abstol,'alpha',alpha,'nmax',nmax)
%
% [pHat, out_param] = *meanMCBer_g*(Yrand,in_param)
%% Description
%
% pHat = *meanMCBer_g*(Yrand) estimates the mean of a Bernoulli random
%  variable Y to within a specified absolute error tolerance with
%  guaranteed confidence level 99%. Input Yrand is a function handle that
%  accepts a positive integer input n and returns a n x 1 vector of IID
%  instances of the Bernoulli random variable Y.
% 
% pHat = *meanMCBer_g*(Yrand,abstol,alpha,nmax) estimates the mean
%  of a Bernoulli random variable Y to within a specified absolute error
%  tolerance with guaranteed confidence level 1-alpha using all ordered
%  parsing inputs abstol, alpha and nmax.
% 
% pHat = *meanMCBer_g*(Yrand,'abstol',abstol,'alpha',alpha,'nmax',nmax)
%  estimates the mean of a Bernoulli random variable Y to within a
%  specified absolute error tolerance with guaranteed confidence level
%  1-alpha. All the field-value pairs are optional and can be supplied in
%  different order.
% 
% [pHat, out_param] = *meanMCBer_g*(Yrand,in_param) estimates the mean
%  of a Bernoulli random variable Y to within a specified absolute error
%  tolerance with the given parameters in_param and produce the estimated
%  mean pHat and output parameters out_param.
% 
% *Input Arguments*
%
% * Yrand --- the function for generating IID instances of a Bernoulli
%            random variable Y whose mean we want to estimate.
%
% * pHat --- the estimated mean of Y.
%
% * in_param.abstol --- the absolute error tolerance, the default value is 1e-2.
% 
% * in_param.alpha --- the uncertainty, the default value is 1%.
% 
% * in_param.nmax --- the sample budget, the default value is 1e9.
% 
% *Output Arguments*
%
% * out_param.n --- the total sample used.
%
% * out_param.time --- the time elapsed in seconds.
% 
%%  Guarantee
%
% If the sample size is calculated according Hoeffding's inequality, which
% equals to ceil(log(2/out_param.alpha)/(2*out_param.abstol^2)), then the
% following inequality must be satisfied:
%
% Pr(|p-pHat| <= abstol) >= 1-alpha.
% 
% Here p is the true mean of Yrand, and pHat is the output of MEANMCBER_G
%
% Also, the cost is deterministic.
%
%%   Examples
%   *Example 1*

%   Calculate the mean of a Bernoulli random variable with true p=1/90,
%   absolute error tolerance 1e-3 and uncertainty 0.01.
% 
    in_param.abstol=1e-3; in_param.alpha = 0.01; p=1/9;Yrand=@(n) rand(n,1)<p;
    pHat = meanMCBer_g(Yrand,in_param)
 
%% 
%   *Example 2*

%   Using the same function as example 1, with the absolute error tolerance
%   1e-4.
% 
    pHat = meanMCBer_g(Yrand,1e-4)
    
%% 
%   *Example 3*

%   Using the same function as example 1, with the absolute error
%   tolerance 1e-2 and uncertainty 0.05.
% 
    pHat = meanMCBer_g(Yrand,'abstol',1e-2,'alpha',0.05)
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
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
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
% Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014.
% arXiv:1208.4318 [math.ST]
%
% [2] Sou-Cheng T.  Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% "GAIL: Guaranteed Automatic Integration Library (Version 2.1)"
% [MATLAB Software], 2015. Available from
% http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software", Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, DOI:
% http://dx.doi.org/10.5334/jors.bb, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Sou-Cheng T. Choi, "Summary of the First Workshop On Sustainable
% Software for Science: Practice And Experiences (WSSSPE1)", Journal of
% Open Research Software, Volume 2, Number 1, e6, pp. 1-21, DOI:
% http://dx.doi.org/10.5334/jors.an, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
