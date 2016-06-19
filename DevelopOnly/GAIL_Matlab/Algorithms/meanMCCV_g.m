function [tmu,out_param]=meanMCCV_g(varargin)
% meanMCCV_G Monte Carlo method to estimate the mean of a random variable.
% with the option to add control variates
%
%   Input Arguments
%
%     YXrand --- the function for generating n IID instances of a pair of
%     random variables (Y, X_c), where
%     Y is the random variable whose mean we want to estimate. X_c is the
%     control variate of dimension p (corresponding to the number of
%     control variates we waent to add). Y is often defined as a function
%     of some random variable X with a simple distribution. X_d should
%     generally be a subset of X, with known means. The input of YXrand
%     should be the number of random variables n, the output of YXrand
%     should be a n-by-(p+1) matrix, with the first column containing the
%     values of Y and the remaining p+1 columns containing the values of
%     the p control variates.
%     For example, ...
%
%`    muX --- a row vector containing the means of all control variates
%
%     in_param.abstol --- the absolute error tolerance, which should be
%     positive, default value is 1e-2.
%
%     in_param.reltol --- the relative error tolerance, which should be
%     between 0 and 1, default value is 1e-1.
%
%     in_param.alpha --- the uncertainty, which should be a small positive
%     percentage. default value is 1%.
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
%     out_para.cv --- control variate status
% 
%                   0   Success
%                   1   Unable to use control variates
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

start = tic; %start the clock
[YXrand, muX, Yrand, randindex, out_param] = meanMCCV_g_param(varargin{:});

Yrandop = Yrand_CV(YXrand, muX, Yrand, randindex, out_param);

nb = out_param.nb*length(muX);

[tmu,out_param]=meanMC_g(Yrandop, out_param);

out_param.nb = nb;
out_param.ntot = out_param.ntot + nb + randindex * 4e3; % update ntot
out_param.time = toc(start); % update time

end

function Yrandop = Yrand_CV(YXrand, muX, Yrand, randindex, out_param)
% this is the control variate part, follows the formula of control variate
n = out_param.nb*length(muX);
YX = YXrand(n);% generate the matrix YX
Y =YX(:,1); % get Y
X =bsxfun(@minus,YX(:,2:end),mean(YX(:,2:end))); % get centered X

b = regress(Y, [ones(n,1) X]); % the optimal beta estimated

    function Y = Yrand_b(n) % new estimator with control variates
        YX_b = YXrand(n);
        Y_b =YX_b(:,1);
        X_b =YX_b(:,2:end);
        Y = Y_b - bsxfun(@minus,X_b,muX)*b(2:end);
    end
Yrandop = @Yrand_b;
if randindex % if randY is provided
    fn1 = @() Yrand(1e3); 
    fn2 = @() Yrand_b(1e3);

    t1 = timeit(fn1); % time to run Yrand 1000 times
    t2 = timeit(fn2); % time to run Yrand_b 1000 times

    var1 = var(Yrand(1e3)); % estimated var of the plain estimator
    var2 = var(Yrand_b(1e3)); % estimated var of the estimator with control variates

    inflation = 1.2; 
    if var2*t2 < var1*t1*inflation % the modified estimator performs better
        Yrandop = @Yrand_b;
    else % the plain estimator performs better
        Yrandop = Yrand;
    end
end
end

function  [YXrand, muX, Yrand, randindex, out_param] = meanMCCV_g_param(varargin)

default.abstol  = 1e-2;% default absolute error tolerance
default.reltol = 1e-1;% default relative error tolerance
default.nb = 1e4; % default initial sample size nSig for estimation of 
                  % control covariate coefficients (per covariate)
default.nSig = 1e4;% default initial sample size nSig for variance estimation
default.n1 = 1e4; % default initial sample size n1 for mean estimation
default.alpha = 0.01;% default uncertainty
default.fudge = 1.2;% default fudge factor
default.tbudget = 100;% default time budget
default.nbudget = 1e9; % default sample budget

randindex = 0; % when Yrand is not provided

if isempty(varargin)
%    help meanMCCV_g
%    warning('GAIL:meanMCCV_g:yrandnotgiven',...
%        'YXrand must be specified. Now GAIL is using YXrand =@(n) rand(n,1).^2.')
    YXrand = @(n) rand(n,1).^2;
    muX = [];
    %if no values are parsed, print warning message and use the default
    %random variable
else
    YXrand = varargin{1};
end

if numel(varargin)<2
%    help meanMCCV_g
%    warning(...)
    muX = [];
else
    muX = varargin{2};
end

index = 3;
if numel(varargin)>=index
    if isa(varargin{index}, 'function_handle')
        Yrand = varargin{3};
        index = index + 1;
        randindex = 1;
    else 
      Yrand = @(n) rand(n,1).^2; 
    end
else 
    Yrand = @(n) rand(n,1).^2; 
end

validvarargin=numel(varargin)>=index;
if validvarargin
    in4=varargin{index};
    validvarargin=(isnumeric(in4) || isstruct(in4) ...
        || ischar(in4));
end

if ~validvarargin
    %if there is only input which is YXrand, use all the default parameters
    out_param.abstol = default.abstol;% default absolute error tolerance
    out_param.reltol = default.reltol; % default relative error tolerance
    out_param.alpha = default.alpha;% default uncertainty
    out_param.fudge = default.fudge;% default standard deviation inflation factor
    out_param.nb = default.nb; % default sample size to estimate control covariate coef
    out_param.nSig = default.nSig;% default sample size to estimate the variance
    out_param.n1 = default.n1;% default the initial sample size to estimate the mean
    out_param.tbudget = default.tbudget;% default time budget
    out_param.nbudget = default.nbudget;% default sample budget
else
    p = inputParser;
    addRequired(p,'YXrand'); % need a validation check function here...
    addRequired(p, 'muX'); % need a validation check function here... 
    addRequired(p, 'Yrand'); % need a validation check function here... 
    if isnumeric(in4)
        %if there are multiple inputs with only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'nb',default.nb,@isnumeric);
        addOptional(p,'nSig',default.nSig,@isnumeric);
        addOptional(p,'n1',default.n1,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
    else
        if isstruct(in4) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'nb',default.nb,@isnumeric);
        addParamValue(p,'nSig',default.nSig,@isnumeric);
        addParamValue(p,'n1',default.n1,@isnumeric);
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
    end
    parse(p,YXrand,muX,Yrand,varargin{index:end})
    out_param = p.Results;
end
%if (~gail.isfcn(YXrand))
%    warning('GAIL:meanMC_g:yrandnotfcn',...
%        ['Yrand must be a function handle.'...
%        ' Now GAIL is using default Yrand =@(n) rand(n,1).^2 .'])
    %print warning message
%    YXrand = @(n) rand(n,1).^2;
%end
%if max(size(YXrand(5)))~=5 || min(size(YXrand(5)))~=1
%    % if the input is not a length n vector, print the warning message
%    warning('GAIL:meanMC_g:yrandnotlengthN',...
%        ['Yrand should be a random variable vector of length n, '...
%        'but not an integrand or a matrix.'...
%        ' Now GAIL is using the default YXrand =@(n) rand(n,1).^2.'])
%    YXrand = @(n) rand(n,1).^2;
%end

%pass the signal indicating the parameters have been checked
end