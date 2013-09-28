% Lan Jiang ljiang14@hawk.iit.edu 
% I want to learn how to use repository and how to make GAIL more reliable
% efficient and rubust. 

clear all;close all;clc
yrand = @(n) exp(rand(n,1).^2);
[mu out_param] = meanMC_g(yrand, 1e-3)

% mu =
% 
%     1.4624
% 
% 
% out_param = 
% 
%                   abstol: 1.0000e-03
%                    alpha: 0.0100
%                    fudge: 1.1000
%                  n_sigma: 1000
%                  nbudget: 100000000
%                   npcmax: 1000000
%               timebudget: 100
%                    Yrand: @(n)exp(rand(n,1).^2)
%           n_left_predict: 105948464
%     time_n_sigma_predict: 0.0075
%                     nmax: 99998975
%                      var: 0.2351
%                  kurtmax: 1.1497
%                     n_mu: 2314690
%                       mu: 1.4624
%                        n: 2315690
%                     time: 0.1428
% 

% when meanMC print the error, it does not show the total sample size used.
% need to modify the code. also, need to combine parsing function input and
% parameter checking.
