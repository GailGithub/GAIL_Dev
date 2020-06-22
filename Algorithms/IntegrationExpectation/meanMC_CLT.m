function [sol, out] = meanMC_CLT(varargin)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   sol = MEANMC_CLT(Y,absTol,relTol,alpha,nSig,inflate) estimates the
%   mean, mu, of a random variable to within a specified error tolerance,
%   i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
%   1-alpha, where abstol is the absolute error tolerance.  The default
%   values are abstol=1e-2 and alpha=1%. Input Y is a function handle that
%   accepts a positive integer input n and returns an n x 1 vector of IID
%   instances of the random variable.
%
%   This is a heuristic algorithm based on a Central Limit Theorem
%   approximation.
%
%
%   Input Arguments
%
%     Y --- the function or structure for generating n IID instances of a
%     random variable Y whose mean we want to estimate. Y is often defined
%     as a function of some random variable X with a simple distribution.
%     The input of Yrand should be the number of random variables n, the
%     output of Yrand should be n function values. For example, if Y = X.^2
%     where X is a standard uniform random variable, then one may define
%     Yrand = @(n) rand(n,1).^2.
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
%     Y --- the random generator
%
%     absTol --- the absolute error tolerance
%
%     relTol --- the relative error tolerance
%
%     alpha --- the uncertainty
%
%     mu --- the estimated mean of Y.
%
%     stddev --- sample standard deviation of the random variable
%
%     nSample --- total sample used.
%
%     time --- the time elapsed in seconds.
%
%     errBd --- the error bound.
%
% _Authors: Yueyi Li, Hu Cauw Hung, Fred J. Hickernell_
%
% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval
% [0,1)^2 with absolute tolerance 1e-3 and relative tolerance 0:
% >> [mu,out] = meanMC_CLT(@(n) prod(rand(n,2),2), 0.001);
% >> exact = 1/4;
% >> check = double(abs(exact - mu) < 2e-3)
% check = 1
%
%
% Example 2:
% Estimate the integral f(x)=exp(-x^2) in the interval [0,1] using x as a
% control variate and relative error 1e-3:
% >> f = @(x)[exp(-x.^2), x];
% >> YXn = @(n)f(rand(n,1));
% >> s = struct('Y',YXn,'nY',1,'trueMuCV',1/2);
% >> exact = erf(1)*sqrt(pi)/2;
% >> success = 0; runs = 1000; tol = 1e-3;
% >> for i=1:runs, success = success + double(abs(exact-meanMC_CLT(s,0,tol)) < tol*exact); end
% >> check = success >= 0.99 * runs
%    1
%
%
% Example 3:
% Estimate the Keister's integration in dimension 1 with a=1, 1/sqrt(2) and
% using cos(x) as a control variate:
% >> normsqd = @(x) sum(x.*x,2);
% >> f = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)).* exp((1/2-a^2)*normt);
% >> f1 = @(x,a,d) f(normsqd(x),a,d);
% >> f2 = @(x)[f1(x,1,1),f1(x,1/sqrt(2),1),cos(x)];
% >> YXn = @(n)f2(randn(n,1));
% >> s = struct('Y',YXn,'nY',2,'trueMuCV',1/sqrt(exp(1)));
% >> [hmu,out] = meanMC_CLT(s,0,1e-3);
% >> exact = 1.380388447043143;
% >> check = double(abs(exact-hmu) < max(0,1e-3*abs(exact)))
% check = 1
%
%
% Example 4:
% Estimate the integral with integrand f(x) = x1.^3.*x2.^3.*x3.^3 in the
% interval [0,1]^3 with pure absolute error 1e-3 using x1.*x2.*x3 as a
% control variate:
%
% >> f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).*x(:,2).*x(:,3)];
% >> s = struct('Y',@(n)f(rand(n,3)),'nY',1,'trueMuCV',1/8);
% >> [hmu,out] = meanMC_CLT(s,1e-3,0);
% >> exact = 1/64;
% >> check = double(abs(exact-hmu) < max(1e-3,1e-3*abs(exact)))
% check = 1
%
%
% Example 5:
% Estimate the integrals with integrands f1(x) = x1.^3.*x2.^3.*x3.^3 and
% f2(x)= x1.^2.*x2.^2.*x3.^2-1/27+1/64 in the interval [0,1]^3 using
% x1.*x2.*x3 and x1+x2+x3 as control variates:
%
% >> f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64,  x(:,1).*x(:,2).*x(:,3), x(:,1)+x(:,2)+x(:,3)];
% >> s = struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5]);
% >> [hmu,out] = meanMC_CLT(s,1e-4,1e-3);
% >> exact = 1/64;
% >> check = double(abs(exact-hmu) < max(1e-4,1e-3*abs(exact)))
% check = 1
%
%
%  References
%
%   [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
% _Authors: Yueyi Li, Cu Hauw Hung, Fred J. Hickernell_

tstart = tic; %start the clock

if nargin
   if isa(varargin{1},'gail.cubMCOut')
      out = varargin{1};
   end
end
if ~exist('out','var')
   out = gail.meanYOut(gail.meanYParam(varargin{:}));
end
Yrand = out.Y; %the random number generator
q = out.nY; %the number of target random variable
p = out.CM.nCV; %the number of control variates
val = Yrand(out.nSig); %get samples to estimate variance
if p==0 && q==1
   YY = val(:,1); 
   theta = 1;
else
%if there is control variate, construct a new random variable that has the
%same expected value and smaller variance
   meanVal = mean(val); %the mean of each column
   A = bsxfun(@minus, val, meanVal); %covariance matrix of the samples
   [U, S, V] = svd([A; [ones(1,q) zeros(1,p)] ],0); %use SVD to solve a constrained least square problem
   Sdiag = diag(S); %the vector of the single values
   U2 = U(end,:); %last row of U
   theta = V*(U2'/(U2*U2')./Sdiag); %get the coefficient for control variates
   YY = [val(:,1:q) A(:,q+1:end)] * theta; %get samples of the new random variable 
end

out.stddev = std(YY); %standard deviation of the new samples
sig0up = out.CM.inflate .* out.stddev; %upper bound on the standard deviation
hmu0 = mean(YY); %mean of the samples

nmu = max(1,ceil((-gail.stdnorminv(out.alpha/2)*sig0up ...
   /max(out.err.absTol,out.err.relTol*abs(hmu0))).^2));
   %number of samples needed for the error tolerance
if nmu > out.CM.nMax %don't exceed sample budget
   warning(['The algorithm wants to use nmu = ' int2str(nmu) ...
      ', which is too big. Using ' int2str(out.CM.nMax) ' instead.'])
   nmu = out.CM.nMax; %revise nmu
end

W = Yrand(nmu); %get samples for computing the mean

if p == 0 && q == 1 %no control variates
   YY = W;
else %samples of the new random variable
  if ~isempty(out.CM.trueMuCV)
     W(:,q+1:end) = bsxfun(@minus, W(:,q+1:end), out.CM.trueMuCV); 
     %subtract true mean from control variates
  end
  YY = W * theta; %incorporate the control variates and multiple Y's
end


sol = mean(YY); %estimated mean
out.sol = sol; %record answer in output class

out.nSample = out.nSig+nmu; %total samples required
out.errBd = -gail.stdnorminv(out.alpha/2)*sig0up/sqrt(nmu);
out.theta = theta; %coefficients of the control variates
out.time = toc(tstart); %elapsed time
end
