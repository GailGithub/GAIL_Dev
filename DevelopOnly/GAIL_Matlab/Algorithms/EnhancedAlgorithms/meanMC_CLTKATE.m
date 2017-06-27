function [hmu,mean_out]=meanMC_CLTKATE(varargin)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   tmu = MEANMC_CLT(Y,absTol,relTol,alpha,nSig,inflate) estimates the
%   mean, mu, of a random variable Y to within a specified error tolerance,
%   i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
%   1-alpha, where abstol is the absolute error tolerance.  The default
%   values are abstol=1e-2 and alpha=1%. Input Yrand is a function handle
%   that accepts a positive integer input n and returns an n x 1 vector of
%   IID instances of the random variable Y.
%
%   Input Arguments
%
%     Y --- the function or structure for generating n IID instances of a random
%     variable Y whose mean we want to estimate. Y is often defined as a
%     function of some random variable X with a simple distribution. The
%     input of Yrand should be the number of random variables n, the output
%     of Yrand should be n function values. For example, if Y = X.^2 where X
%     is a standard uniform random variable, then one may define Yrand =
%     @(n) rand(n,1).^2.
%     
%
%     absTol --- the absolute error tolerance, which should be
%     non-negative --- default = 1e-2
%
%     relTol --- the relative error tolerance, which should be
%     non-negative and no greater than 1 --- default = 0
%
%     alpha --- the uncertainty, which should be a small positive
%     percentage --- default = 1%
%
%     nSig --- the number of samples used to compute the sample variance
%     --- default = 1000
%
%     inflate --- the standard deviation inflation factor --- default = 1.2
%
%   Output Arguments
%
%     hmu --- the estimated mean of Y.
%
%     out_param.ntot --- total sample used.
%
%     out_param.var --- the sample variance.
%
%     out_param.time --- the time elapsed in seconds.
%
% >> [mu,out] = meanMC_CLT(@(n) rand(n,1).^2, 0.001)
% mu =
%     0.33***
% out = 
%   meanYOut with properties:
% 
%           mu: 0.33***
%          std: 0.***
%            Y: @(n)rand(n,1).^2
%        alpha: 0.0100
%         nSig: 1000
%      inflate: 1.2000
%         nMax: 100000000
%       absTol: 1.0000e-03
%       relTol: 0
%       solFun: @(mu)mu
%     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
%      nSample: ***
%         time: ***
%
%



% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2 with absolute 
% tolerance 1e-5 and relative tolerence 0:
% 
% >> f = @(x) prod(x,2);
% >> q = meanMC_CLTKATE(@(n)f(rand(n,2)),1e-5,0); exactsol = 1/4; 
% >> check = abs(exactsol-q) < 1e-5
% check = 1
%
%



% Example 2:
% Estimate the integral with integrand f(x) = x1.^3.*x2.^3.*x3.^3
% in the interval [0,1)^3 with pure absolute error 1e-3 using x1.*x2.*x3 as control variate:
% 
% >> f=@(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).*x(:,2).*x(:,3)];
% >> s=struct('Y',@(n)f(rand(n,3)),'nY',1,'trueMuCV',1/8)
% >> [hmu,mean_out]=meanMC_CLTKATE(s,1e-3,0) exactsol = 1/64;
% >> check = abs(exactsol-hmu) < max(1e-3,1e-3*abs(exactsol))
% check = 1
%
% Example 3:
% Estimate the integrals with integrands f1(x) = x1.^3.*x2.^3.*x3.^3 and 
% f2(x)= x1.^2.*x2.^2.*x3.^2-1/27+1/64 in the interval [0,1)^3
% using x1.*x2.*x3 and x1+x2.^3+x3 as control variate:
%
% >> f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64,x(:,1).*x(:,2).*x(:,3),x(:,1)+x(:,2)+x(:,3)];
% >> s=struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5])
% >> [hmu,mean_out]=meanMC_CLTKATE(s,1e-4,1e-3); exactsol = 1/64;
% >> check = abs(exactsol-hmu) < max(1e-4,1e-3*abs(exactsol))
% check = 1
%
% Example 4:
% Estimate the Keister's integration in dimension 1 with a=1, 1/sqrt(2)and using cos(x) as a control variate:
% >> normsqd = @(t) sum(t.*t,2);
% >> f=@(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)).* exp((1/2-a^2)*normt);
% >> f1 = @(t,a,d) f(normsqd(t),a,d);
% >> f2=@(t)[f1(t,1,1),f1(t,1/sqrt(2),1),cos(t)];
% >> YXn=@(n)f2(randn(n,1));
% >> s=struct('Y',YXn,'nY',2,'trueMuCV',1/sqrt(exp(1)))
% >> [hmu,mean_out]=meanMC_CLTKATE(s,0,1e-3); exactsol=1.380388447043143;
% >>  abs(exactsol-hmu) < max(0,1e-3*abs(exactsol))




% This is a heuristic algorithm based on a Central Limit Theorem
% approximation

tstart = tic; %start the clock 
mean_inp = gail.meanYParam(varargin{:}); %parse the input and check it for errors
mean_out = gail.meanYOut(mean_inp); %create the output class
Yrand=mean_out.Y; %the random number generator
q=mean_out.nY; %the number of target random variable 
p=mean_out.nCV; %the number of control variates
xmean=mean_out.trueMuCV; %the mean of the control variates

val = Yrand(mean_out.nSig); %get samples to estimate variance 
if p==0 && q==1
    YY = val(:,1); 
    
else
%if there is control variate, construct a new random variable that has the
%same expected value and smaller variance
        meanVal=mean(val); %the mean of each column
        A=bsxfun(@minus, val, meanVal); %covariance matrix of the samples
        [U,S,V]=svd([A; [ones(1,q) zeros(1,p)] ],0); %use SVD to solve a constrained least square problem
        Sdiag = diag(S); %the vector of the single values
        U2=U(end,:); %last row of U
        beta=V*(U2'/(U2*U2')./Sdiag); %get the coefficient for control variates
        YY = [val(:,1:q) A(:,q+1:end)] * beta; %get samples of the new random variable 
end

mean_out.std = std(YY); %standard deviation of the new samples

sig0up = mean_out.inflate .* mean_out.std; %upper bound on the standard deviation
hmu0 = mean(YY); % mean of the samples

nmu = max(1,ceil((-gail.stdnorminv(mean_out.alpha/2)*sig0up ...
   /max(mean_out.absTol,mean_out.relTol*abs(hmu0))).^2)); 
   %number of samples needed for the error tolerance
if nmu > mean_out.nMax %don't exceed sample budget
   warning(['The algorithm wants to use nmu = ' int2str(nmu) ...
      ', which is too big. Using ' int2str(mean_out.nMax) ' instead.']) 
   nmu = mean_out.nMax;
end

YY = Yrand(nmu); %get samples 

if p > 0 || q > 1   %samples of the new random variable
  YY(:,q+1:end) = bsxfun(@minus, YY(:,q+1:end), xmean);
  YY = YY*beta;
end

hmu = mean(YY); %estimated mean
mean_out.mu = hmu;
mean_out.nSample = mean_out.nSig+nmu; %total samples required
mean_out.time = toc(tstart); %elapsed time
end


