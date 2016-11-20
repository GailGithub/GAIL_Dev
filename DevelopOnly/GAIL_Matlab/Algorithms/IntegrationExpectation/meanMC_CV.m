function [tmu,out_param]=meanMC_CV(varargin)
% YOPTPRICE_CV creates the control variate output for option pricing using
% the |optPayoff| object
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

if isempty(varargin)
    help meanMC_CV
    warning('GAIL:meanMC_CV:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand =@(n) [rand(n,1).^2, rand(n,1)].')
    YXrand = @(n) [rand(n,1).^2, rand(n,1)];
    %if no values are parsed, print warning message and use the default
    %random variable
else
    YXrand = varargin{1};
end

if numel(varargin)>1
    meanX = varargin{2};
else
    help meanMC_CV
    warning('GAIL:meanMC_CV:meanXnotgiven',...
        'meanX must be given. Now GAIL is using sample mean of X.')
    YX = YXrand(10000);
    X = YX(:,2:end);
    meanX = mean(X,1);
end

beta =@(YXrand) (bsxfun(@minus,YXrand(:,2:end),mean(YXrand(:,2:end),1))\YXrand(:,1)); %optimal beta
YCVrand =@(YXrand) (YXrand(:,1) - bsxfun(@minus,YXrand(:,2:end),meanX)*beta(YXrand)); %control variate random variable
[tmu, out_param] = meanMC_g(@(n) (YCVrand(YXrand(n))), varargin{3:end});