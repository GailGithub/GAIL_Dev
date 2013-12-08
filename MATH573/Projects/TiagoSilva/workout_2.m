% In this second workout, we want to estimate the convergence rate of the
% IPaS estimator by analyzing the relationship between standard deviation 
% of the estimation and the sample size M. Since we don't know how to 
% estimate the standard deviation of the estimation only by using the same 
% sample used to estimate the true value, then we have to estimate de 
% standard deviation empirically by running the estimator 100 times. Even 
% though workout1 showed that the IPaS estimator might be more efficient 
% than vanilla MC method, this workout shows that both has the same 
% convergence rate, i.e. log(StDev) is proportional to (1/2)*log(M).
%
% Here, we fixed the IPaS levels at scenario.split=[1,2,3,4,5,6,7];
% using 10 bernoulli trials to estimate the probobability that the sum of 
% successful outcomes  is greather than 7. For each sample size M, we
% calculate the StDev of the estimation.

clear all;close all;clc;

N = 8;             % Number of scenarios
Ns = N*(N+1)/2+N;   % Total number of simulations

scenario.split=[1,2,3,4,5,6,7]; scenario.M=2000;
gamma = FindExactSolForBinoProblem(scenario.split(end),10,0.1);

x=zeros(Ns,1);
y=zeros(Ns,1);
k=0;
for i=1:N
    for j=1:(i+1)
    k=k+1;
    x(k)=scenario.M;
    y(k)= EstimateStd(scenario)/gamma;
    end
    scenario.M=ceil(scenario.M / 1.41);
    disp([num2str(100*i/N) '% Done...']);
end

logY = log(y');
logX = log(x');
C = [ones(length(logX'),1) logX']\logY';
alpha = C(2);
error = C(1)+C(2)*logX-logY;
SME = sum(error.^2)/(Ns-2);
SS = sum((logX-mean(logX)).^2);
sb = sqrt(SME/SS);
CI = tinv(0.975,Ns-2)*sb;
disp(['The 95% confidence interval of the slope is  [' num2str(alpha-CI) ',' num2str(alpha+CI) ']'])
hold off
loglog(x,y,'.')
