% This is the driver script to test the meanMC_g algorithm
clear all; close all; clc;
%y = @(n) rand(n,1).^2;% the test function
y = @GAIL_Internal.Ytrafficmodel; % this is the traffic model
in_param.abstol = 1e-2;% the absolute error tolerance
% in_param.alpha = 0.01;% uncertainty
% in_param.n_sigma = 1e3;% sample size to estimate the variance
% in_param.fudge = 1.1;% variance inflation factor
% in_param.timebudget = 10;% time budget
in_param.nbudget = 1e8/3;% sample budget
% in_param.npcmax = 1e6;% optimal piesewise maximum to calculate mu
[mu, out_param] = meanMC_g(y,in_param) % the results

%% The following output was obtain on 2014-6-25
% % 
% 
% mu =
% 
%    4.459808444599791
% 
% 
% out_param = 
% 
%                   abstol: 0.010000000000000
%                    alpha: 0.010000000000000
%                  checked: 2
%                    fudge: 1.100000000000000
%                  n_sigma: 1000
%                  nbudget: 100000000
%                   npcmax: 1000000
%                  tbudget: 100
%                    Yrand: @GAIL_Internal.Ytrafficmodel
%           n_left_predict: 13363
%     time_n_sigma_predict: 8.224439100000000
%                     nmax: 13363
%                      var: 0.008407935937481
%                  kurtmax: 1.149741494296920
%                     n_mu: 2861
%                       mu: 4.459808444599791
%                        n: 3887
%                     time: 33.299483418999998
% 
% 
