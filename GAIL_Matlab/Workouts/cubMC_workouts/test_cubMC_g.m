% This is the drive script to test cubMC_g algorithm
clear all;close all;clc;
in_param.measure  = 'uniform';
%f=@(x) exp(-x(1).^2-x(2).^2);% the test function
%f=@(x) x(:,1)+x(:,2);
in_param.dim =3;%the function dimension
startingpoint = zeros(1,in_param.dim);
endingpoint = ones(1,in_param.dim);
interval = [startingpoint;endingpoint];% the integration interval
in_param.abstol = 1e-3;% the absolute tolerance
in_param.alpha = 0.01;% the uncertainty
in_param.n_sigma = 1e4;% the sample size to estimate sigma
in_param.fudge =1.1;% standard deviation inflation factor
in_param.timebudget = 100;% time budget
in_param.nbudget = 1e9;% sample budget
%[Q,out_param]=cubMC_g(f,interval,'alpha',0.05)% the results
index = 2;% the dimension of the test function
alpha = ones(1:in_param.dim); % one coefficient
beta = [1/3 1/4 2];
%beta = ones(1:in_param.dim);% the other coefficent
r=1;
test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
f_true = genz_test_fun_true (interval,index,in_param.dim,alpha,beta)
[Q,out_param]=cubMC_g(test_function,interval,in_param)% the results
