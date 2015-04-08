function [tmu,out_param]=meanMC_CLT(Yrand,abstol,alpha)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   tmu = MEANMC_CLT(Yrand) estimates the mean, mu, of a random variable Y to
%   within a specified error tolerance, i.e., | mu - tmu | <= abstol with
%   probability at least 1-alpha, where abstol is the absolute error
%   tolerance.  The default values are abstol=1e-2 and alpha=1%. Input
%   Yrand is a function handle that accepts a positive integer input n and
%   returns an n x 1 vector of IID instances of the random variable Y.
%
%   Input Arguments
%
%     Yrand --- the function for generating n IID instances of a random
%     variable Y whose mean we want to estimate. Y is often defined as a
%     function of some random variable X with a simple distribution. The
%     input of Yrand should be the number of random variables n, the output
%     of Yrand should be n function values. For example, if Y = X.^2 where X
%     is a standard uniform random variable, then one may define Yrand =
%     @(n) rand(n,1).^2.
%
%     abstol --- the absolute error tolerance, which should be
%     positive, default value is 1e-2.
%
%     alpha --- the uncertainty, which should be a small positive
%     percentage. default value is 1%.
%
%   Output Arguments
%
%     tmu --- the estimated mean of Y.
%
%     out_param.ntot --- total sample used.
%
%     out_param.var --- the sample variance.
%
%     out_param.time --- the time elapsed in seconds.
%

%This is a heuristic algorithm based on a Central Limit Theorem
%approximation
if nargin < 3
   alpha = 0.01;
   if nargin < 2
      abstol = 0.01;
      if nargin < 1
         Yrand = @(n) rand(n,1);
      end
   end
end
nSig=1e4;
nMax=1e8;
out_param.alpha = alpha;
out_param.fudge = 1.2;
tstart = tic; %start the clock
Yval = Yrand(nSig);% get samples to estimate variance 
out_param.var = var(Yval);% calculate the sample variance--stage 1
sig0 = sqrt(out_param.var);% standard deviation
sig0up = out_param.fudge.*sig0;% upper bound on the standard deviation
alpha1 = 1-sqrt(1-out_param.alpha);% the uncertainty for variance estimation
nmu = max(1,ceil((-norminv(alpha1)*sig0up/abstol).^2));
assert(nmu<nMax,['nmu = ' int2str(nmu) ', which is too big'])
tmu=mean(Yrand(nmu));
out_param.ntot=nSig+nmu;
out_param.time=toc(tstart); %elapsed time
end

