% This is the driver script to test the meanMC_g algorithm
clear all; close all; clc;
%y = @(n) rand(n,1).^2;% the test function
y = @Ytrafficmodel; % this is the traffic model
in_param.abstol = 1e-2;% the absolute error tolerance
in_param.reltol = 0;
in_param.tbudget = 50;
in_param.nsig = 1e2;
in_param.n1 = 1e2;
in_param.alpha = 0.05;% uncertainty
in_param.fudge = 1.1;
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
