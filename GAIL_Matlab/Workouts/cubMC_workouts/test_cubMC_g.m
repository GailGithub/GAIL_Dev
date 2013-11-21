% This is the drive script to test cubMC_g algorithm
clear all;close all;clc;
in_param.measure  = 'uniform';
%f=@(x) exp(-x(1).^2-x(2).^2);% the test function
%f=@(x) x(:,1)+x(:,2);
in_param.dim = 2;%the function dimension
interval = [zeros(1,in_param.dim);ones(1,in_param.dim)];% the integration interval
in_param.abstol = 1e-2;% the absolute tolerance
in_param.alpha = 0.01;% the uncertainty
in_param.n_sigma = 1e4;% the sample size to estimate sigma
in_param.fudge =1.1;% standard deviation inflation factor
in_param.timebudget = 100;% time budget
in_param.nbudget = 1e9;% sample budget
%[Q,out_param]=cubMC_g(f,interval,'alpha',0.05)% the results
index = 2;% the dimension of the test function
alpha = [1 2]; % one coefficient
beta = [1 2 ];% the other coefficent
r=1;
test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
[Q,out_param]=cubMC_g(test_function,interval,'alpha',0.05)% the results
