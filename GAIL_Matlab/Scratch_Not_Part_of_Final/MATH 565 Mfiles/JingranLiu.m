%% Test 2 Jingran's method
clear all
%reset(RandStream.getDefaultStream);
n=10;

%% Control variates
n=1000 %sample size
f=@(x) exp(x(:,1).*x(:,2)); %define function to be integrated
%f=@(x) exp(x(:,1)+x(:,2)); %define function to be integrated
x=rand(n,2); %uniform random samples
cv=x(:,1).*x(:,2).^2;
mucv=1/6; %true means of x_1 and x_2
cvbar=mean(cv,1); %sample averages of x_1 and x_2
G=cv-cvbar; %x - xbar in matrix form
y=f(x); %function values
ybar=mean(y) %sample mean of function values
stdy=std(y)
ciMCwidth=1.96*std(y)/sqrt(n) %confidence interval for simple MC
beta=G\(y-ybar); %beta vector for control variates
muhat=ybar+(mucv-cvbar)*beta %control variate estimator
resid=y-ybar-G*beta; %residuals from regression
stddev=std(resid) %standard error of residuals
ciCVwidth=1.96*stddev/sqrt(n) %confidence interval width for control variates
errratio=ciMCwidth/ciCVwidth %ratio of errors

