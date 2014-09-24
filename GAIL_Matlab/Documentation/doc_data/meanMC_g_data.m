%  Guarantee
% This algorithm attempts to calculate the mean of a random variable to a
% prescribed error tolerance with guaranteed confidence level 1-alpha. If
% the algorithm terminated without showing any warning messages and provide
% an answer tmu, then the follow inequality would be satisfied:
% 
% Pr(|mu-tmu| <= max(abstol,reltol|mu|)) >= 1-alpha
%
% where abstol is the absolute error tolerance and reltol is the relative
% error tolerance, if the true mean mu is rather small as well as the
% reltol, then the abstol would be satisfied, and vice versa. 
%
% The cost of the algorithm is also bounded above by N_up, which is defined
% in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And the
% following inequality holds:
% 
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
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
