function [tmu,out_param]=meanMC_CV(varargin)
% meanMC_CV Monte Carlo method estimate the mean of a random variable by
% using control variate.
%
%   tmu = meanMC_CV(YXrand) estimates the mean, mu, of a random variable Y 
%   by inputting other depentent random variables X to within a specified 
%   generalized error tolerance.
%   tolfun:=max(abstol,reltol*|mu|), i.e., |mu - tmu| <= tolfun with
%   probability at least 1-alpha, where abstol is the absolute error
%   tolerance, and reltol is the relative error tolerance. Usually the
%   reltol determines the accuracy of the estimation, however, if the |mu|
%   is rather small, the abstol determines the accuracy of the estimation.
%   The default values are abstol=1e-2, reltol=1e-1, and alpha=1%. Input
%   YXrand is a function handle that accepts a positive integer input n and
%   returns an n x j array of IID instances of the random variable Y and
%   the depentent random variables X.
%
%   tmu = meanMC_CV(YXrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%   estimates the mean of a random variable Y by inputting other depentent
%   random variables X to within a specified generalized error tolerance 
%   tolfun with guaranteed confidence level 1-alpha using all ordered parsing 
%   inputs abstol, reltol, alpha, fudge, nSig, n1, tbudget, nbudget.
%
%   tmu = meanMC_CV(YXrand,'abstol',abstol,'reltol',reltol,'alpha',alpha,
%   'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',tbudget,'nbudget',nbudget)
%   estimates the mean of a random variable Y by inputting other depentent
%   random variables X to within a specified generalized error tolerance 
%   tolfun with guaranteed confidence level 1-alpha. All the field-value pairs 
%   are optional and can be supplied in different order, if a field is not 
%   supplied, the default value is used.
%
%   [tmu, out_param] = meanMC_CV(YXrand,in_param) estimates the mean of a
%   random variable Y by inputting other depentent random variables X to 
%   within a specified generalized error tolerance tolfun with the given 
%   parameters in_param and produce the estimated mean tmu and output parameters 
%   out_param. If a field is not specified, the default value is used.
%
%   Input Arguments
%
%     YXrand --- the function for generating n IID instances of a random
%     variable Y whose mean we want to estimate and other dependent random
%     variables X. Y is often defined as a function of some random variable
%     with a simple distribution. The input of Yrand should be the number 
%     of random variables n, the output of Yrand should be n x j function 
%     values. For example, if Y = u.^2 and X = u where u is a standard 
%     uniform random variable, then one may define 
%     YXrand = @(n) [rand(n,1).^2, rand(n,1)].
%
%     MeanX --- the population mean of X, default value is the sample mean
%     of X.
%
%     in_param.abstol --- the absolute error tolerance, which should be
%     positive, default value is 1e-2.
%
%     in_param.reltol --- the relative error tolerance, which should be
%     between 0 and 1, default value is 1e-1.
%
%     in_param.alpha --- the uncertainty, which should be a small positive
%     percentage, default value is 1%.
%
%     in_param.fudge --- standard deviation inflation factor, which should
%     be larger than 1, default value is 1.2.
%
%     in_param.nSig --- initial sample size for estimating the sample
%     variance, which should be a moderate large integer at least 30, the
%     default value is 1e4.
%
%     in_param.n1 --- initial sample size for estimating the sample mean,
%     which should be a moderate large positive integer at least 30, the
%     default value is 1e4.
%
%     in_param.tbudget --- the time budget in seconds to do the two-stage
%     estimation, which should be positive, the default value is 100 seconds.
%
%     in_param.nbudget --- the sample budget to do the two-stage
%     estimation, which should be a large positive integer, the default
%     value is 1e9.
%
%   Output Arguments
%
%     tmu --- the estimated mean of Y.
%
%     out_param.tau --- the iteration step.
%
%     out_param.n --- the sample size used in each iteration.
%
%     out_param.nremain --- the remaining sample budget to estimate mu. It was
%     calculated by the sample left and time left.
%
%     out_param.ntot --- total sample used.
%
%     out_param.hmu --- estimated mean in each iteration.
%
%     out_param.tol --- the reliable upper bound on error for each iteration.
%
%     out_param.var --- the sample variance.
%
%     out_param.exit --- the state of program when exiting.
%
%                      0   Success
%
%                      1   Not enough samples to estimate the mean
%
%     out_param.kurtmax --- the upper bound on modified kurtosis.
%
%     out_param.time --- the time elapsed in seconds.
%
%     out_param.flag --- parameter checking status
%
%                           1  checked by meanMC_g
%
%  Guarantee
% This algorithm attempts to calculate the mean, mu, of a random variable
% to a prescribed error tolerance, tolfun:= max(abstol,reltol*|mu|), with
% guaranteed confidence level 1-alpha. If the algorithm terminated without
% showing any warning messages and provide an answer tmu, then the follow
% inequality would be satisfied:
%
% Pr(|mu-tmu| <= tolfun) >= 1-alpha
%
% The cost of the algorithm, N_tot, is also bounded above by N_up, which is
% defined in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And
% the following inequality holds:
%
% Pr (N_tot <= N_up) >= 1-beta
%
%
% Examples

if isempty(varargin)
    % if no values are parsed, print warning message and use the default
    % random variable
    help meanMC_CV
    warning('GAIL:meanMC_CV:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand =@(n) [rand(n,1).^2, rand(n,1)].')
    YXrand = @(n) [rand(n,1).^2, rand(n,1)];
else
    YXrand = varargin{1};
end

if numel(varargin)>1 && isnumeric(varargin{2})
    meanX = varargin{2};
else
    % if no meanX values are parsed, print warning message and using sample
    % mean of X
    help meanMC_CV
    warning('GAIL:meanMC_CV:meanXnotgiven',...
        'meanX must be given. Now GAIL is using sample mean of X.')
    YX = YXrand(10000);
    X = YX(:,2:end);
    meanX = mean(X,1);
end

% optimal beta
beta =@(YXrand) (bsxfun(@minus,YXrand(:,2:end),mean(YXrand(:,2:end),1))\...
    YXrand(:,1));
% control variate random variable
YCVrand =@(YXrand) (YXrand(:,1) - bsxfun(@minus,YXrand(:,2:end),meanX)...
    *beta(YXrand)); 
% use meanMC_g to compute mu
[tmu, out_param] = meanMC_g(@(n) (YCVrand(YXrand(n))), varargin{3:end});