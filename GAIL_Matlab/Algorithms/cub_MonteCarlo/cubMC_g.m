function [Q,out_param] = cubMC_g(f,interval,varargin)
% CUBMC_G Monte Carlo method to evaluate a multidimentional integral to
% within a specified absolute error tolerance with guaranteed uncertainy
% within alpha.
%
%   [Q,out_param] = CUBMC_G(f) estimates the integral with integrand f to
%   within the absolute error tolerance 1e-2 and with guaranteed
%   uncertainty alpha within 1%. Input f a function handle. The function
%   Y=f(X) should accept a vector argument X and return a vector result Y,
%   the integrand evaluated at each element of X.
%
%   Q = CUBMC_G(f,interval,measure,abstol,alpha,n_sigma,fudge) estimates the
%   integral with integrand f to within an absolute error tolerance abstol
%   with guaranteed uncertainty within alpha using ordered parameter input
%   interval, measure, tolerence, uncertainty, n_sigma and fudge factor.
%
%   Q =
%   CUBMC_G(f,interval,'measure','uniform','abstol',abstol,'alpha',alpha,
%   'n_sigma',n_sigma,fudge',fudge) estimates the integral with integrand f
%   to within an absolute error tolerance abstol with guaranteed
%   uncertainty within alpha. All the field-value pairs are optional and
%   can be supplied in different order.
%
%   Q = CUBMC_G(f,interval,in_param) estimates the integral with integrand f
%   to within an absolute error tolerance in_param.abstol with guaranteed
%   uncertainty within in_param.alpha. If a field is not specified, the
%   default value is used.
%
%   f --- the integrand.
%
%   interval --- the integration interval.
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
%   out_param.time --- the time eclipsed.
%
%   out_param.exit --- the state of program when exiting.
%                         0   success
%                         10  interval does not contain numbers
%                         11  interval not 2 x d
%                         12  interval is only a point in one direction
%                         13  interval is infinite when measure is uniform
%                         14  interval is not doubly infinite when measure is normal
%   Examples:
%
%   Example 1:
%   Estimate the integral with integrand f(x) = x.^2 in the interval [1,2]
%   
%   >> f=@(x) x.^2;interval = [1;2];
%   >> Q = cubMC_g(f,interval,'abstol',1e-2)
%   Q = 2.33***
%
%
%   Example 2:
%   Estimate the integral with integrand f(x) = x.^2 in the interval [1,2]
%
%   >> f=@(x) x.^2;interval = [1;2];
%   >> Q = cubMC_g(f,interval)
%   Q = 2.33***
%
%
%   See also FUNAPPX_G, INTEGRAL_G, MEANMC_G
%
%   Reference:
%   [1]  Hickernell, F. J., Jiang, L., Liu, Y. and Owen, A. B., 
%        Guaranteed Conservative Fixed Width Confidence Intervals 
%        Via Monte Carlo Sampling, preprint, 2013, 	arXiv:1208.4318 [math.ST]. 

tstart=tic;
default.measure = 'uniform';% default measure
default.dim = 1;% default dimension
default.interval = [zeros(1,default.dim);ones(1,default.dim)];% default interval
default.abstol  = 1e-2;% default absolute error tolerence
default.alpha = 0.01;% default uncertainty
default.n_sigma = 1e3; % default n_sigma
default.fudge = 1.1; % default variance inflation factor

if (nargin<1) % if no input print error message
    help cubMC_g
    warning('f must be specified. Now GAIL is using f = @(x) x.^2.')
    f = @(x) x.^2;
end;

if (nargin<2) 
% if only one input which is function f, use default values for all other
% parameters
    interval = default.interval;
    in_param.measure = default.measure;
    in_param.abstol = default.abstol;
    in_param.alpha = default.alpha;
    in_param.n_sigma = default.n_sigma;
    in_param.fudge = default.fudge;
end;
p = inputParser;
addRequired(p,'f',@isfcn);

if (nargin<3) 
% if two inputs which are functon f and interval, then use all the default
% parameters for all other inputs.
    in_param.measure = default.measure;
    in_param.abstol = default.abstol;
    in_param.alpha = default.alpha;
    in_param.n_sigma = default.n_sigma;
    in_param.fudge = default.fudge;
end;
p = inputParser;
addRequired(p,'f',@isfcn);
addRequired(p, 'interval', @isnumeric);

if (nargin == 3 && isstruct(varargin{1})) 
    % add f interval and in_param as input
    p.StructExpand = true;
    p.KeepUnmatched  = true;
    addParamValue(p,'measure',default.measure,@(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
    addParamValue(p,'abstol',default.abstol,@isnumeric);
    addParamValue(p,'alpha',default.alpha,@isnumeric);
    addParamValue(p,'n_sigma',default.n_sigma,@isnumeric);
    addParamValue(p,'fudge',default.fudge,@isnumeric);   
    parse(p,f,interval,varargin{:})
    in_param = p.Results;
end

if (nargin > 3)
    in2 = varargin{1}; % the input parameters are not in order
    if (ischar(in2)),
        addParamValue(p,'measure',default.measure,@(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'n_sigma',default.n_sigma,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);        
        parse(p,f,interval,varargin{:})
        in_param = p.Results;
    end
end

if (nargin >= 3 && isnumeric(varargin{1})) % the input parameters are in order
    addOptional(p,'measure',default.measure,@(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
    addOptional(p,'abstol',default.abstol,@isnumeric);
    addOptional(p,'alpha',default.alpha,@isnumeric);
    addOptional(p,'n_sigma',default.n_sigma,@isnumeric);
    addOptional(p,'fudge',default.fudge,@isnumeric);    
    parse(p,f,interval,varargin{:})
    in_param = p.Results;
end
in_param = cubMC_g_param(interval,default,in_param);%check validity of inputs
out_param = in_param; %let the out_param contains all the in_param
f=transformIntegrand(f,interval,in_param); 
%transform integrand so that the interval would not need to be changed
if strcmp(in_param.measure,'uniform')%the using uniformly distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(rand(nfun,in_param.dim)),in_param);
    out_param.Q=Q;% using meanMC_g to get the mean 
else strcmp(in_param.measure,'normal')%using normally distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(randn(nfun,in_param.dim)),in_param);
    out_param.Q=Q; % using meanMC_g to get the mean
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

function [in_param,out_param] = cubMC_g_param(interval,default,in_param)

out_param.exit=0; %success! until found otherwise

[two, in_param.dim]=size(interval); %interval should be 2 x dimension
if two==0 && isfield(in_param,'interval'); %if interval specified through in_param structure
    interval=in_param.interval; %then get it from there
    [two, in_param.dim]=size(interval); %and get the dimension
end
if any(isnan(interval(:))); %check interval for not a number
    out_param.exit=10; out_param = cubMC_g_err(out_param); return; 
end
if two~=2 %if interval is given as row vector for dimension 1, fix that
    if in_param.dim==2; in_param.dim=two; interval=interval'; 
    else out_param.exit=11; out_param = cubMC_g_err(out_param); return; %else return an error
    end
end
interval=[min(interval,[],1); max(interval,[],1)]; %ensure left and right endpoints are in order
if any(interval(1,:)==interval(2,:)); %interval is a point in one direction
    out_param.exit=12; out_param = cubMC_g_err(out_param); return;
end
in_param.interval=interval; %copy interval into the param structure

if isfield(in_param,'measure');
    in_param.measure=validatestring(in_param.measure,{'uniform','normal','Gaussian'});
    if strcmp(in_param.measure,'Gaussian'); in_param.measure='normal'; end
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
if (~isposint(in_param.n_sigma)) %Initial sample size
    warning('the number n_sigma should a positive integer,')
    warning('take the absolute value and ceil.')
    in_param.n_sigma = ceil(abs(in_param.n_sigma));
end
if (in_param.fudge <= 1) %Variance inflation factor/fudge factor
    warning('the fudge factor should be bigger than 1, ')
    warning('use the default value.')
    in_param.fudge = default.fudge;
end
if (in_param.abstol <= 0) %error tolerence
    warning('the absolute error tolerence should be larger than 0, ')
    warning('use the absolute value.')
    in_param.abstol = abs(in_param.abstol);
end
if (in_param.alpha <= 0 ||in_param.alpha >= 1) %Uncertainty 
    warning('the uncertainy should be less than 1 and bigger than 0, ')
    warning('use the the default value.')
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
    case 10; fprintf(2,'Error: interval must contain numbers\n');
    case 11; fprintf(2,'Error: interval must be 2 x d\n');
    case 12; fprintf(2,...
        'Error: interval must be more than a point in any coordinate direction\n');
    case 13; fprintf(2,'Error: interval must be finite when measure is uniform\n');
    case 14; fprintf(2,['Error: interval must be infinite in both directions' ...
        ' when measure is normal\n']);
end
out_param.Q=NaN;
Q=out_param.Q;
if nargin>1; out_param.time=toc(tstart); end
end