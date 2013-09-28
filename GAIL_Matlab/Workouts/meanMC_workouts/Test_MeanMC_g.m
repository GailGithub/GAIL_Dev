% This is the drive file to test the meanMC_g algorithm
clear all; close all; clc;
y = @(n) rand(n,1).^2;% the test function
%y = @Ytrafficmodel; % this is the traffic model
in_param.abstol = 1e-2;% the absolute error tolerance
in_param.alpha = 0.01;% uncertainty
in_param.n_sigma = 1e3;% sample size to estimate the variance
in_param.fudge = 1.1;% variance inflation factor
in_param.timebudget = 100;% time budget
in_param.nbudget = 1e8;% sample budget
in_param.npcmax = 1e6;% optimal piesewise maximum to calculate mu
[mu, out_param] = meanMC_g(y,in_param) % the results
