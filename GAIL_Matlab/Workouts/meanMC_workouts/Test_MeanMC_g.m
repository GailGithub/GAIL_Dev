% This is the driver script to test the meanMC_g algorithm
clear all; close all; clc;
%y = @(n) rand(n,1).^2;% the test function
y = @Ytrafficmodel; % this is the traffic model
in_param.abstol = 1e-2;% the absolute error tolerance
% in_param.alpha = 0.01;% uncertainty
% in_param.n_sigma = 1e3;% sample size to estimate the variance
% in_param.fudge = 1.1;% variance inflation factor
% in_param.timebudget = 10;% time budget
% in_param.nbudget = 1e8;% sample budget
% in_param.npcmax = 1e6;% optimal piesewise maximum to calculate mu
[mu, out_param] = meanMC_g(y,in_param) % the results

%% The following output was obtain on 2014-2-10
% 
% mu =
% 
%     4.4557
% 
% 
% out_param = 
% 
%                   abstol: 0.0100
%                    alpha: 0.0100
%                  checked: 2
%                    fudge: 1.1000
%                  n_sigma: 1000
%                  nbudget: 100000000
%                   npcmax: 1000000
%                  tbudget: 100
%                    Yrand: @Ytrafficmodel
%           n_left_predict: 13063
%     time_n_sigma_predict: 8.0745
%                     nmax: 13063
%                      var: 0.0089
%                  kurtmax: 1.1497
%                     n_mu: 2988
%                       mu: 4.4557
%                        n: 4014
%                     time: 37.0792
