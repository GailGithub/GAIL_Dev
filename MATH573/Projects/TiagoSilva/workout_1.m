% In this first workout, we want to compare the efficiency of three
% differente estimators. First estimator is based on vanilla Monte Carlo
% method. The second estimator is based on the IPaS method, using 3
% split levels, and the third estimator is again based on IPaS method, but
% using 6 split levels. In this analysis, we fix the tolerance os the
% estimation by the following: The standard deviation of the estimation is
% 50% of the estimation. By doing that, we can compare the sample size
% needed and the computational time necessary to satisfy the same
% condition. That is a good way to analyse performance.
%
% For each case, we want to estimate the probability gamma that X_T, a sum 
% of T bernoulli trials, is greather or equal to 7, where standard 
% deviation of the estimation is at least 50% of its value gamma. Each 
% bernoulli trial has a probability of success equal to 0.1

clear all;close all;clc;

case1.T=10;                   % Default - Number of steps
case1.f = @(v,U) v+(U<0.1);   % Default - Markov function f:X_i -> X_(i+1)
case2=case1; case3=case1;     % Default - case2 and case3  

case1.split = [7];           % estimator 1 - naive MC
case2.split = [2,4,7];       % estimator 2 - using IPaS with 3 levels
case3.split = [2,3,4,5,6,7]; % estimator 3 - using IPaS with 6 levels

tol = 0.5;                    % tol = std(gamma_hat)/gamma
gamma = FindExactSolForBinoProblem(7,10,0.1); % Exact solution for gamma
case1.M=ceil((1/tol)^2*(1-gamma)/gamma); % Sample size that respect tol
sig = sqrt(gamma*(1-gamma))/sqrt(case1.M); % Exact variance of the estimator 1


%% Estimate gamma using Estimator 1
[output1.gamma,output1.elapsed_time]=IPaS(case1); output1.M=case1.M;
output1.stdev = sqrt(output1.gamma*(1-output1.gamma))/sqrt(output1.M);
disp('20% Done...');
%% Find M2 and M3, for estimator 2 and estimator 3 respectively, that 
% deliver aproximately the same variance for the estimation when comparing 
% to estimator 1
x2=zeros(100,1); x3=zeros(100,1);
M_test=5*10^2;              % Sample Size test

std2_test=EstimateStd(case2.f,case2.split,case2.T,M_test);
disp('40% Done...');
std3_test=EstimateStd(case3.f,case3.split,case3.T,M_test);
disp('60% Done...');
%Here we have the same size case2.M and case3.M, for Estimator 2 and 
%Estimator 3 respectively, that delivers and estimation with variance 
%similar to our Estimator 1
case2.M = ceil(M_test *(std2_test/sig)^2);
case3.M = ceil(M_test *(std3_test/sig)^2);

%% Using case2.M and case3.M find estimator of gamma using estimator 2 and 3
et2=zeros(100,1); et3=zeros(100,1);
[output2.stdev,output2.gamma, output2.elapsed_time]= EstimateStd(case2);
disp('80% Done...');
[output3.stdev,output3.gamma, output3.elapsed_time]= EstimateStd(case3);
disp('100% Done...');

output2.M=case2.M; 
output3.M=case3.M; 


%% Output
disp('For these three simple examples:')
disp('case1 - lvl=7; case2 - lvl=[2,4,7]; case3 - lvl=[2,3,4,5,6,7]')
disp(['True Value = ' num2str(gamma)])
disp(['Estimation case1 = ' num2str(output1.gamma)])
disp(['Estimation case2 = ' num2str(output2.gamma)])
disp(['Estimation case3 = ' num2str(output3.gamma)])
disp(['Standard Dev case1 = ' num2str(output1.stdev)])
disp(['Standard Dev case2 = ' num2str(output2.stdev)])
disp(['Standard Dev case3 = ' num2str(output3.stdev)])
disp(['Sample size case1 = ' num2str(output1.M)])
disp(['Sample size case2 = ' num2str(output2.M)])
disp(['Sample size case3 = ' num2str(output3.M)])
disp(['Average Time case1 = ' num2str(output1.elapsed_time)])
disp(['Average Time case2 = ' num2str(output2.elapsed_time)])
disp(['Average Time case3 = ' num2str(output3.elapsed_time)])


