%% Release Notes 
%
% GAIL Version 2.3.1, 2020.
% 
%% Major changes in algorithms
%
% In this release, we have a new algorithm called CubBayesNet_g. Similar to
% cubBayesLattice_g, it is an automatic Bayesian cubature that considers
% the integrand a realization of a Gaussian process. CubBayesNet_g uses
% Sobol points whereas cubBayesLattice_g uses lattice points.
% 
% 
%% Major changes in publications
% 
% In the folder "Papers", we have added a few recently published research
% articles and a thesis related to our core algorithms.
%
% First, we have Rathinavel's 2019 PhD thesis and his joint publication
% with Hickernell that develop the theory behind cubBayesLattice_g and
% CubBayesNet_g.
%
% In addition, we have included Ding, Hickernell, and Jimenez Rugama's
% recent paper, An Adaptive Algorithm Employing Continuous Linear
% Functionals.
% 
% 
%% Bug Fixes
% 
% In the previous versions of GAIL, three of our multiple integration
% algorithms, cubMC_g, cubLattice_g, and cubSobol_g, when applied to an
% integral on a non-unit hypercube with respect to uniform measure, the
% numerical approximations omitted the proper normalization constant, i.e.,
% dividing by the volume of the hypercube. This problem has been resolved
% in this release. We illustrate the bug fix with a simple example.
% Consider f(x,y) = exp(-x^2-y^2) with (x,y) in D:=[-1,2]^2. Integrating
% f with respect to the Lebesgue measure, we obtain I ~= 2.65333.
% Integrating f with respect to the uniform measure, we have instead I / 9,
% since the integration domain has volume 9. Our GAIL functions by default
% integrate with respect to the uniform measure, but previous versions
% returned answers with respect to the Lebesgue measure.