function [Q,out_param] = cubMC_g(varargin)
% CUBMC_G Monte Carlo method to evaluate a multidimensional integral to
% within a specified absolute error tolerance with guaranteed uncertainty
% within alpha. 
%
%  The guarantee holds if the modified kurtosis is less than the kmax,
%  which is defined in terms of uncertainty(alpha), sample size to estimate
%  variance(n_sigma) and standard deviation inflation factor(fudge). For
%  details, please refer to our paper.
%
%   [Q,out_param] = CUBMC_G(f,interval) estimates the integral with
%   integrand f to within the absolute error tolerance 1e-2 and with
%   guaranteed uncertainty alpha within 1%. Input f a function handle. The
%   function Y=f(X) should accept a vector argument X and return a vector
%   result Y, the integrand evaluated at each element of X. Input interval
%   is 2 x d matrix. 
%
%   Q = CUBMC_G(f,interval,measure,abstol,alpha,n_sigma,fudge) estimates the
%   integral with integrand f to within an absolute error tolerance abstol
%   with guaranteed uncertainty within alpha using ordered parameter input
%   interval, measure, tolerance, uncertainty, n_sigma and fudge factor.If
%   an input is not specified, the default value is used.
%
%   Q =
%   CUBMC_G(f,interval,'measure','uniform','abstol',abstol,'alpha',alpha,
%   'n_sigma',n_sigma,fudge',fudge) estimates the integral with integrand f
%   to within an absolute error tolerance abstol with guaranteed
%   uncertainty within alpha. All the field-value pairs are optional and
%   can be supplied in different order. If an input is not specified, the
%   default value is used.
%
%   Q = CUBMC_G(f,interval,in_param) estimates the integral with integrand f
%   to within an absolute error tolerance in_param.abstol with guaranteed
%   uncertainty within in_param.alpha. If a field is not specified, the
%   default value is used.
%
%   f --- the integrand.
%
%   interval --- the integration interval. The default interval is [0 1]^d.
%
%   in_param.measure --- the measure for generating the random variable,
%   the default is uniform.
%
%   in_param.abstol --- the absolute error tolerance, default value is 1e-2.
%
%   in_param.alpha --- the uncertainty, default value is 1%.
%
%   in_param.n_sigma --- initial sample size for estimating the sample
%   variance, the default value is 1e3.
%
%   in_param.fudge --- the standard deviation inflation factor, the default
%   value is 1.1.
%   
%   Q --- the estimated value of the the integration.
%
%   out_param_time_n_sigma_predict --- the estimated time to get n_sigma
%   samples of the random variable.
%
%   out_param.n_left_predict --- using the time left to predict the number
%   of samples left.
%
%   out_param.nmax --- the maximum sample budget to estimate mu, it comes
%   from both the sample budget and the time budget.
%
%   out_param.var --- the sample variance.
%
%   out_param.kurtmax --- the upper bound on modified kurtosis.
%
%   out_param.time --- the time eclipsed.
%
%   out_param.n_mu --- the sample size that needed to estimate the mu.
%
%   out_param.n --- the total sample size needed to do the two stage
%   algorithm.
%
%   out_param.exit --- the state of program when exiting.
%                         0   success
%                         1   No enough samples to estimate the mean.
%                         2   Initial try out time costs more than
%                             10% of time budget. 
%                         3   The estimated time for estimating variance 
%                             is bigger than half of the time budget.
%                         10  Interval does not contain numbers.
%                         11  Interval not 2 x d.
%                         12  Interval is only a point in one direction.
%                         13  Interval is infinite when measure is uniform.
%                         14  Interval is not doubly infinite when measure
%                             is normal.
%   Examples:
%
%   Example 1:
%   Estimate the integral with integrand f(x) = x.^2 in the interval [0,1]
%   
%   >> f=@(x) x.^2;interval = [0;1];
%   >> Q = cubMC_g(f,interval,'abstol',1e-3)
%   Q = 0.33***
%
%
%   Example 2:
%   Estimate the integral with integrand f(x) = exp(x) in the interval [1,2]
%
%   >> f=@(x) exp(x);interval = [1;2];
%   >> Q = cubMC_g(f,interval,'uniform',1e-3)
%   Q = 4.67***
%
%
%   Example 3:
%   Estimate the integral with integrand f(x) = sin(x) in the interval [1,2]
%
%   >> f=@(x) sin(x);interval = [1;2];
%   >> Q = cubMC_g(f,interval,'uniform',1e-3)
%   Q = 0.95***
%
%
%   Example 4: 
%   Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
%   interval [0 0;1 1],where x is a vector x = [x1 x2].
%
%   >> f=@(x) exp(-x(:,1).^2-x(:,2).^2);interval = [0 0;1 1];
%   >> Q = cubMC_g(f,interval,'uniform',1e-3)
%   Q = 0.55***
%
%
%   Example 5: 
%   Estimate the integral with integrand f(x) = 2^n*prod(x1*x2*...*xn)+0.555 in the
%   interval [zeros(1,n);ones(1,n)],where x is a vector x = [x1 x2 ... xn].
%
%   >> n=3;f=@(x) 2^n*prod(x,2)+0.555;interval = [zeros(1,n);ones(1,n)];
%   >> Q = cubMC_g(f,interval,'uniform',1e-3)
%   Q = 1.5***
%
%
%   See also FUNAPPX_G, INTEGRAL_G, MEANMC_G
%
%   Reference:
%   [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G.
%   W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to
%   appear, arXiv:1208.4318 [math.ST]

tstart=tic;
[f,interval,in_param,out_param] = cubMC_g_param(varargin{:});%check validity of inputs
f=transformIntegrand(f,interval,in_param); 
% transform integrand so that the interval would not need to be changed
if strcmp(in_param.measure,'uniform')% the using uniformly distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(rand(nfun,in_param.dim)),in_param);
    out_param.Q=Q;% using meanMC_g to get the mean 
else strcmp(in_param.measure,'normal')% using normally distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(randn(nfun,in_param.dim)),in_param);
    out_param.Q=Q;% using meanMC_g to get the mean
end
out_param.time=toc(tstart); %elapsed time
end
function newf=transformIntegrand(oldf,interval,in_param) 
    if strcmp(in_param.measure,'uniform') %uniform measure
        a=interval(1,:); %left endpoint
        b=interval(2,:); %right endpoint
    if all(a==0) && all(b==1) %no change needed
        newf=oldf; 
    else %transform points and integrand
        bmina=b-a; %interval width
        volbox=prod(bmina); %volume of the interval
        newf=@(x) oldf(x.*repmat(bmina,size(x,1),1)+repmat(a,size(x,1),1))...
            .*volbox;
       %stretch and shift, then multiply by volume
    end
    elseif strcmp(in_param.measure,'normal')
        newf=oldf;% no change if it is normal measure.
    end   
end

function [f,interval,in_param,out_param] = cubMC_g_param(varargin)

default.measure = 'uniform';% default measure
default.dim = 1;% default dimension
default.interval = [zeros(1,default.dim);ones(1,default.dim)];% default interval
default.abstol  = 1e-2;% default absolute error tolerence
default.alpha = 0.01;% default uncertainty
default.n_sigma = 1e3; % default n_sigma
default.fudge = 1.1; % default variance inflation factor

if isempty(varargin) % if no input print error message and use the default setting
    help cubMC_g
    warning('MATLAB:cubMC_g:fnotgiven',['f must be specified. Now GAIL is using f = @(x) x.^2. '...
        'Integration interval must be specified. Now GAIL is using interval [0 1]'])
    f = @(x) x.^2;
    interval = default.interval;
elseif numel(varargin)==1
    % if there is only function but no interval input. Use default interval.
    help cubMC_g
    warning('MATLAB:cubMC_g:intervalnotgiven',...
        'the interval must be specified, Now GAIL is using interval [0 1]')
    f = varargin{1};
    interval = default.interval;    
else
    f = varargin{1};
    interval = varargin{2}; % the first input is function, the second input is interval.
end
    
validvarargin=numel(varargin)>2;% check if there is any optional parameter input
if validvarargin
    in3=varargin{3}; % check the third input
    validvarargin=(isnumeric(in3) || isstruct(in3) || ischar(in3));
    % to see if it is numeric structure or character.
end

if ~validvarargin
% if there is no optional input, use default settings.
    in_param.measure = default.measure;
    in_param.abstol = default.abstol;
    in_param.alpha = default.alpha;
    in_param.n_sigma = default.n_sigma;
    in_param.fudge = default.fudge;
else % if there is some optional input 
    p = inputParser;
    addRequired(p,'f',@isfcn);
    addRequired(p,'interval',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        %if there are multiple inputs with only numeric, they should be put
        %in order.
        addOptional(p,'measure',default.measure,@(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'n_sigma',default.n_sigma,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
    else
        if isstruct(in3) %the input is structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end % if there are multiple inputs with name and numeric, they
    % could be put not in order
        addParamValue(p,'measure',default.measure,@(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'n_sigma',default.n_sigma,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
    end
    parse(p,f,interval,varargin{3:end})
    in_param = p.Results;
end
out_param = in_param; % let the out_param contains all the in_param
[two, in_param.dim]=size(interval); %interval should be 2 x dimension
if two==0 && isfield(in_param,'interval'); 
    %if interval specified through in_param structure
    interval=in_param.interval; %then get it from there
    [two, in_param.dim]=size(interval); %and get the dimension
end
if any(isnan(interval(:))); %check interval for not a number
    out_param.exit=10; out_param = cubMC_g_err(out_param); return; 
end
if two~=2 %if interval is given as row vector for dimension 1, fix that
    if in_param.dim==2; in_param.dim=two; interval=interval';
    else out_param.exit=11; out_param = cubMC_g_err(out_param); return; 
        %else, return an error
    end
end
interval=[min(interval,[],1); max(interval,[],1)]; 
%ensure left and right endpoints are in order
if any(interval(1,:)==interval(2,:)); %interval is a point in one direction
    out_param.exit=12; out_param = cubMC_g_err(out_param); return;
end
in_param.interval=interval; %copy interval into the param structure

if isfield(in_param,'measure'); % the sample measure
    in_param.measure=validatestring(in_param.measure,{'uniform','normal','Gaussian'});
    if strcmpi(in_param.measure,'Gaussian')
        in_param.measure='normal'; 
    end
else
    in_param.measure=default.measure;
end
if strcmp(in_param.measure,'uniform')&&~all(isfinite(interval(:)))
    %cannot integrate on an infinite interval with the uniform distribution
    out_param.exit=13; out_param = cubMC_g_err(out_param); return;
end
if strcmp(in_param.measure,'normal')&&any(isfinite(interval(:)))
    %must integrate on an infinite interval with the normal distribution
    out_param.exit=14; out_param = cubMC_g_err(out_param); return;
end
if (~isposint(in_param.n_sigma)) %the sample to estimate sigma
    warning('MATLAB:cubMC_g:nsignotposint',...
        ['the number n_sigma should a positive integer,'...
        'take the absolute value and ceil.'])
    in_param.n_sigma = ceil(abs(in_param.n_sigma));
end
if (in_param.fudge <= 1) %standard deviation inflation factor/fudge factor
    warning('MATLAB:cubMC_g:fudgelessthan1',...
        'the fudge factor should be bigger than 1, use the default value.')
    in_param.fudge = default.fudge;
end
if (in_param.abstol <= 0) %error tolerance
    warning('MATLAB:cubMC_g:abstolneg',...
        'the absolute error tolerence should be larger than 0, use the absolute value.')
    in_param.abstol = abs(in_param.abstol);
end
if (in_param.alpha <= 0 ||in_param.alpha >= 1) %uncertainty 
    warning('MATLAB:cubMC_g:alphanot01',...
        ['the uncertainy should be less than 1 and bigger than 0, '...
    'use the the default value.'])
    in_param.alpha = default.alpha;
end

end
function [out_param,Q]=cubMC_g_err(out_param,tstart)
%Handles errors in cubMC_g and cubMC_g_param
%to give an exit with information
%out_param.exit = 0   success
%             10  interval does not contain numbers
%             11  interval not 2 x d
%             12  interval is only a point in one direction
%             13  interval is infinite when measure is uniform
%             14  interval is not doubly infinite when measure is normal
if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
    case 10; error('MATLAB:cubMC_g:intervalnotnum',...
            'interval must contain numbers.');
    case 11; error('MATLAB:cubMC_g:intervalnot2d',...
            'interval must be 2 x d.');
    case 12; error('MATLAB:cubMC_g:intervalnotlessthan2',...
            'interval must be more than a point in any coordinate direction.');
    case 13; error('MATLAB:cubMC_g:intervalnotfiniteforuniform',...
            'interval must be finite when measure is uniform.');
    case 14; error('MATLAB:cubMC_g:intervalnotinffornormal',...
            ['interval must be infinite in both directions' ...
        ' when measure is normal']);
end
out_param.Q=NaN;
Q=out_param.Q;
if nargin>1; out_param.time=toc(tstart); end
end
