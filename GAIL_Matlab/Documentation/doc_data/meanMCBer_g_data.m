%%  Guarantee
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
%%   Examples
%   *Example 1*

%   Calculate the mean of a Bernoulli random variable with true p=1/90,
%   absolute error tolerance 1e-3 and uncertainty 0.01.
% 
    in_param.abstol=1e-3; in_param.alpha = 0.01; p=1/9;Yrand=@(n) rand(n,1)<p;
    pHat = meanMCBernoulli_g(Yrand,in_param)
 
%% 
%   *Example 2*

%   Using the same function as example 1, with the relative error tolerance
%   1e-2.
% 
    pHat = meanMCBernoulli_g(Yrand,0,1e-2,'rel')
    
%% 
%   *Example 3*

%   Using the same function as example 1, with the relative error
%   tolerance 1e-2 and uncertainty 0.05.
% 
    pHat = meanMCBernoulli_g(Yrand,'errtype','rel','reltol',1e-2,'alpha',0.05)
