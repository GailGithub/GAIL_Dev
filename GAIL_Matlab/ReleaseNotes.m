% Release Notes
% 
% GAIL Version 2.3, 2019.
% 
% (1) Major changes in algorithms
% 
% First, we have a new algorithm called cubBayesLattice_g, an efficient
% automatic Bayesian cubature that considers the integrand a realization of
% a Gaussian process.
% 
% In addition, integral_g for univariate integration is enhanced by
% adopting the more efficient Simpson's rule in place of the trapezoidal
% rule. We deprecated integral_g in version 2.2 and renamed it in this
% version as integralTrap_g.
% 
% 
% (2) Major changes in publications
% 
% In the folder "Papers", we have added a few recently published research
% articles and theses related to our core algorithms.
% 
% First of all, it is the doctoral thesis by Yizhi Zhang that develops the
% theory behind integral_g.
% 
% Next, we have Jagadeeswaran and Hickernell (2018+), which elaborates on
% the theory behind the development of cubBayesLattice_g.
% 
% Lastly, a short paper by Hickernell et al. (2018) reviews and compares
% the three GAIL algorithms for high-dimensional integration based on
% (quasi-)Monte Carlo methods, namely cubMC_g, cubLattice_g, and
% cubSobol_g.
% 
% 
% (3) Major changes in tests
% 
% We continue to execute automated nightly fast tests and weekly long tests
% on our server. Moreover, these tests are now conducted for all MATLAB
% versions from R2016a to R2019a. The test reports are available on Mega
% cloud storage at https://mega.nz/. More specifically, fast and long test
% reports are archived at
%     https://mega.nz/#F!4olmWa6L!vYuscSnGqkvkZrGJXW5Umw
% and
%     https://mega.nz/#F!I0cAEKJD!AyQ_8tmxkknfIsuEW0_jnA
% respectively.
% 
