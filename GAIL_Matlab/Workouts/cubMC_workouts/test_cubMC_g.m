% This is the drive script to test cubMC_g algorithm
clear all;close all;clc;
in_param.measure  = 'uniform';
f=@(x) exp(-x.^2*1e6);% the test function
interval = [-0.001;0.001];% the integration interval
in_param.abstol = 1e-2;% the absolute tolerance
in_param.alpha = 0.01;% the uncertainty
in_param.n_sigma = 1e4;% the sample size to estimate sigma
in_param.fudge =1.1;% standard deviation inflation factor
in_param.timebudget = 100;% time budget
in_param.nbudget = 1e9;% sample budget
[Q,out_param]=cubMC_g(f,interval,in_param)% the results
%[Q,out_param] = cubMC_g(f,interval,'abstol',1e-3)