% This is the driver script to test the meanMC_g algorithm
%clear all; close all; clc;
%y = @(n) rand(n,1).^2;% the test function
function [mu,out_param] = Test_meanMC_g
y = @Ytrafficmodel; % this is the traffic model
%y = @(n) rand(n,1).^5;% the test function
in_param.abstol = 1e-2;% the absolute error tolerance
in_param.reltol = 1e-1;%the relative error tolerance
in_param.tbudget = 50;%the time budget
in_param.nSig = 1e2;%the sample size to estimate the variance
in_param.n1 = 1e2;%the initial sample size to estimate the mean
in_param.alpha = 0.01;% uncertainty
in_param.fudge = 1.2;%standard deviation inflation factor

[mu, out_param] = meanMC_g(y,in_param); % the results

end

%% The following output was obtain on 2014-9-30
% 
% mu =
% 
%     4.4604
% 
% 
% out_param = 
% 
%      abstol: 0.0100
%       alpha: 0.0500
%       fudge: 1.1000
%          n1: 100
%     nbudget: 1.0000e+09
%        nSig: 100
%      reltol: 0
%     tbudget: 50
%       Yrand: @Ytrafficmodel
%     checked: 1
%        exit: 0
%        nmax: 6173
%         var: 0.0107
%     kurtmax: 1.0570
%         tol: [0.0931 0.0080]
%           n: [100 4094]
%         tau: 2
%         hmu: [4.4502 4.4604]
%        ntot: 4316
%        time: 32.9001
