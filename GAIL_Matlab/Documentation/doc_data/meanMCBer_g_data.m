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
%%
% *Example 1*

% Calculate the mean of a Bernoulli random variable with true p=1/90,
% absolute error tolerance 1e-3 and uncertainty 0.01.
 
    in_param.abstol=1e-3; in_param.alpha = 0.01;
    p=1/9; Yrand=@(n) rand(n,1)<p;
    pHat = meanMCBer_g(Yrand,in_param)
 
%%
% *Example 2*

% Using the same function as example 1, with the absolute error tolerance
% 1e-4.

    pHat = meanMCBer_g(Yrand,1e-4)
    
%%
% *Example 3*

% Using the same function as example 1, with the absolute error tolerance
% 1e-2 and uncertainty 0.05.

    pHat = meanMCBer_g(Yrand,'abstol',1e-2,'alpha',0.05)
