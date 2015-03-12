% This is the driver script to test the meanMC_g algorithm
clear all; close all; clc;
%y = @(n) rand(n,1).^2;% the test function
%function [mumeanMC,mumeanMCabs,mumeanMCnew] = Test_meanMCnew_g
y = @Ytrafficmodel; % this is the traffic model
%y = @(n) rand(n,1).^2/10;% the test function
in_param.abstol = 1e-2;% the absolute error tolerance
in_param.reltol =0;%the relative error tolerance
in_param.tbudget = 50;%the time budget
in_param.nSig = 1e4;%the sample size to estimate the variance
in_param.n1 = 1e4;%the initial sample size to estimate the mean
in_param.alpha = 0.01;% uncertainty
in_param.nbudget = 1e10;
in_param.fudge = 1.2;%standard deviation inflation factor
tic
[mumeanMC, out_param_meanMC] = meanMC_g(y,in_param); % the results
t_meanMC_g = toc;
tic
[mumeanMCabs, out_param_meanMCabs] = meanMCabs_g(y,in_param); % the results
t_meanMCabs_g = toc;
tic
[mumeanMCnew, out_param_meanMCnew] = meanMCnew_g(y,in_param); % the results
t_meanMCnew_g = toc;
tic
[mumeanMCnew2, out_param_meanMCnew2] = meanMCnew2_g(y,in_param); % the results
t_meanMCnew2_g = toc;
timevec = [t_meanMC_g t_meanMCabs_g t_meanMCnew_g t_meanMCnew2_g]

%end
