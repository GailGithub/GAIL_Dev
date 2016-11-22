function [tmu,out_param]=meanMC_CV(varargin)
% meanMC_CV Monte Carlo method estimate the mean of a random variable through 
% control variate.
%
%   tmu = meanMC_CV(YXrand,meanX) estimates the mean, mu, of a random variable Y 
%   by adding other depentent random variables X to within a specified 
%   generalized error tolerance.
%   tolfun:=max(abstol,reltol*|mu|), i.e., |mu - tmu| <= tolfun with
%   probability at least 1-alpha, where abstol is the absolute error
%   tolerance, and reltol is the relative error tolerance. Usually the
%   reltol determines the accuracy of the estimation, however, if the |mu|
%   is rather small, the abstol determines the accuracy of the estimation.
%   The default values are abstol=1e-2, reltol=1e-1, and alpha=1%. Input
%   YXrand is a function handle that accepts a positive integer input n and
%   returns an n x j matrix of IID instances of the random variable Y and
%   the depentent random variables X.
%
%   tmu = meanMC_CV(YXrand,meanX,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%   estimates the mean of a random variable Y by adding other depentent
%   random variables X to within a specified generalized error tolerance 
%   tolfun with guaranteed confidence level 1-alpha using all ordered parsing 
%   inputs abstol, reltol, alpha, fudge, nSig, n1, tbudget, nbudget.
%
%   tmu = meanMC_CV(YXrand,meanX,'abstol',abstol,'reltol',reltol,'alpha',alpha,
%   'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',tbudget,'nbudget',nbudget)
%   estimates the mean of a random variable Y by adding other depentent
%   random variables X to within a specified generalized error tolerance 
%   tolfun with guaranteed confidence level 1-alpha. All the field-value pairs 
%   are optional and can be supplied in different order, if a field is not 
%   supplied, the default value is used.
%
%   [tmu, out_param] = meanMC_CV(YXrand,meanX,in_param) estimates the mean of a
%   random variable Y by adding other depentent random variables X to 
%   within a specified generalized error tolerance tolfun with the given 
%   parameters in_param and produce the estimated mean tmu and output parameters 
%   out_param. If a field is not specified, the default value is used.
%
%   Input Arguments
%
%     YXrand --- the function for generating n IID instances of a random
%     variable Y whose mean we want to estimate and other dependent random
%     variables X. Y is often defined as a function of some random variable
%     with a simple distribution. The input of YXrand should be the number 
%     of random variables n, the output of YXrand should be n x j function 
%     values. For example, if Y = u.^2 and X = u where u is a standard 
%     uniform random variable, then one may define 
%     YXrand = @(n) [rand(n,1).^2, rand(n,1)].
%
%     MeanX --- the population mean of X, which must be given, 
%     default value is the sample mean of X. If YXrand have been inputted  
%     and the population mean of X is not known, this parameter must be 
%     inputted as not numbers to ensure the result is correct.
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
%
% Example 1:
% If no parameters are parsed, help text will show up and the program will 
% compute the default YXrand: YXrand =@(n) [rand(n,1).^2, rand(n,1)].
%
% >> meanMC_CV
% Warning: YXrand must be specified. Now GAIL is using 
% YXrand =@(n) [rand(n,1).^2, rand(n,1)]. 
% Warning: meanX must be given. Now GAIL is using sample mean of X.
% ans =0.33***
%
%
% Example 2:
% Calculate the mean of u^2 by adding X=u when u is standard normally distributed
% in [-00 +00], with the absolute error tolerance = 1e-3 and uncertainty 5%.
%
% >> in_param.reltol=0; in_param.abstol = 1e-3;
% >> in_param.alpha = 0.05; YXrand=@(n) [randn(n,1).^2, randn(n,1)].
% >> tmu=meanMC_CV(YXrand,0,in_param)
% tmu = 1.00***
%
%
% Example 3:
% Calculate the mean of exp(u) by adding X=u when u is uniformly distributed in
% [0 1], with the absolute error tolerance 1e-3.
%
% >> tmu=meanMC_CV(@(n) [exp(rand(n,1)), rand(n,1)],0.5,1e-3,0)
% tmu = 1.71***
%
%
% Example 4:
% Calculate the mean of cos(u^2) by adding X=[u,u^2] when u is uniformly distributed in
% [0 1], with the relative error tolerance 1e-2 and uncertainty 0.05.
%
% >> tmu=meanMC_CV(@(n) [cos(rand(n,1).^2),rand(n,1), rand(n,1).^2],[0.5,1/3],
% 'reltol',1e-2,'abstol',0,'alpha',0.05)
% tmu = 0.90***

[YXrand,meanX,out_param] = meanMC_CV_param(varargin{:});

% optimal beta
beta =@(YXrand) (bsxfun(@minus,YXrand(:,2:end),mean(YXrand(:,2:end),1))\...
    YXrand(:,1));
% control variate random variable
YCVrand =@(YXrand) (YXrand(:,1) - bsxfun(@minus,YXrand(:,2:end),meanX)...
    *beta(YXrand)); 
% use meanMC_g to compute mu
[tmu, out_param] = meanMC_g(@(n) (YCVrand(YXrand(n))), out_param);
end


function  [YXrand2,meanX2,out_param] = meanMC_CV_param(varargin)

default.abstol  = 1e-2;% default absolute error tolerance
default.reltol = 1e-1;% default relative error tolerance
default.nSig = 1e4;% default initial sample size nSig for variance estimation
default.n1 = 1e4; % default initial sample size n1 for mean estimation
default.alpha = 0.01;% default uncertainty
default.fudge = 1.2;% default fudge factor
default.tbudget = 100;% default time budget
default.nbudget = 1e9; % default sample budget

if isempty(varargin)
    % if no values are parsed, print warning message and use the default
    % random variable
    help meanMC_CV
    warning('GAIL:meanMC_CV:yxrandnotgiven',...
    ['YXrand must be specified. Now GAIL is using the default' ...
    'YXrand =@(n) [rand(n,1).^2, rand(n,1)].'])
    YXrand2 = @(n) [rand(n,1).^2, rand(n,1)];
else
    YXrand2 = varargin{1};
end

if numel(varargin)>1 && isnumeric(varargin{2}) && all(~isnan(varargin{2}))
    YX = YXrand2(100000);
    X = YX(:,2:end);
    meanX2 = mean(X,1);
    if all(abs(meanX2-varargin{2})<0.01) && all(varargin{2}~=1e-1) && ...
            all(varargin{2}~=1e-2) && all(varargin{2}~=1e-3) && ...
            all(varargin{2}~=1e-4) && all(varargin{2}~=1e-5)
        meanX2 = varargin{2};
    else
    % if the difference between any population means and any samples means
    % is not under tolerance 0.01, print warning message and using sample
    % means of X
    warning('GAIL:meanMC_CV:meanXnotright',...
        ['meanX mignt not be given or is not close to sample means of X '...
        'under tolerance 0.01. Now GAIL is using sample means of X.'])
    meanX2 = mean(X,1);
    end
else
    % if meanX are not numbers, print warning message and using sample
    % means of X
    help meanMC_CV
    warning('GAIL:meanMC_CV:meanXnotgiven',...
        'meanX must be given. Now GAIL is using sample means of X.')
    YX = YXrand2(10000);
    X = YX(:,2:end);
    meanX2 = mean(X,1);
end

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin{3};
    validvarargin=(isnumeric(in3) || isstruct(in3) ...
        || ischar(in3));
end

if ~validvarargin
    %if there is only input which is YXrand, use all the default parameters
    out_param.abstol = default.abstol;% default absolute error tolerance
    out_param.reltol = default.reltol; % default relative error tolerance
    out_param.alpha = default.alpha;% default uncertainty
    out_param.fudge = default.fudge;% default standard deviation inflation factor
    out_param.nSig = default.nSig;% default the sample size to estimate the variance
    out_param.n1 = default.n1;% default the initial sample size to estimate the mean
    out_param.tbudget = default.tbudget;% default time budget
    out_param.nbudget = default.nbudget;% default sample budget
else
    p = inputParser;
    addRequired(p,'Yrand',@gail.isfcn);
    if isnumeric(in3)
        %if there are multiple inputs with only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'nSig',default.nSig,@isnumeric);
        addOptional(p,'n1',default.n1,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
    else
        if isstruct(in3) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'nSig',default.nSig,@isnumeric);
        addParamValue(p,'n1',default.n1,@isnumeric);
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
    end
    parse(p,YXrand2,varargin{3:end})
    out_param = p.Results;
end
if (~gail.isfcn(YXrand2))
    %print warning message
    warning('GAIL:meanMC_CV:yxrandnotfcn',...
        ['YXrand must be a function handle.' ...
        ' Now GAIL is using the default YXrand =@(n) [rand(n,1).^2, rand(n,1)].'])
    YXrand2 = @(n) [rand(n,1).^2, rand(n,1)];
end
if size(YXrand2(5),1)~=5 || size(YXrand2(5),2)<1
    % if the input is not a length n x j vector, print the warning message
    warning('GAIL:meanMC_CV:yxrandnotlengthN',...
        ['YXrand should be a random variable n x j matrix, '...
        'but not an integrand.'...
        ' Now GAIL is using the default YXrand =@(n) [rand(n,1).^2, rand(n,1)].'])
    YXrand2 = @(n) [rand(n,1).^2, rand(n,1)];
end
end