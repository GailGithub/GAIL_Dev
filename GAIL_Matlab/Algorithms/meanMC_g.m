function [mu,out_param]=meanMC_g(varargin)
% MEANMC_G Monte Carlo method to estimate the mean of a random variable to
% within a specific absolute error tolerance with guaranteed uncertainty
% within alpha.
%
%   mu = MEANMC_G(Yrand) estimates the mean of a random variable Y to
%   within a specified absolute error tolerance 1e-2 with guaranteed
%   uncertainty within 1%. Input Yrand is a function handle that accepts a
%   positive integer input n and returns an n x 1 vector of IID instances
%   of the random variable Y.
%
%   mu =
%   MEANMC_G(Yrand,abstol,alpha,n_sigma,fudge,timebudget,nbudget,npcmax)
%   estimates the mean of a random variable Y to within an absolute error
%   tolerance abstol with guaranteed uncertainty within alpha using all
%   ordered parsing inputs n_sigma, fudge, timebudget, nbudget and npcmax.
%
%   mu = MEANMC_G(Yrand,'abstol',abstol,'alpha',alpha,'n_sigma',n_sigma,...
%   'fudge',fudge,'timebudget',timebudget,'nbudget',nbudget,'npcmax',npcmax)
%   estimates the mean of a random variable Y to within a specified
%   absolute error tolerance abstol with guaranteed untertainty within
%   alpha. All the field-value pairs are optional and can be supplied in
%   different order.
%
%   mu = MEANMC_G(Yrand,in_param) estimates the mean of a random variable Y
%   to within a specified absolute error tolerance in_param.abstol with
%   guaranteed uncertainty within in_param.alpha. If a field is not
%   specified, the default value is used.
%
%   [mu, out_param] = MEANMC_G(Yrand,in_param) estimates the mean of a
%   random variable Y to within a specified absolute error tolerance with
%   the given parameters in_param and output parameters out_param.
%
%
%   Yrand --- the function for generating IID instances of a random
%   variable Y whose mean we want to estimate. Y is often defined as a
%   function of some random variable X with a simple distribution.  For
%   example, if Y = X.^2 where X is a standard uniform random variable,
%   then one may define Yrand = @(n) rand(n,1).^2.
%
%   mu --- the estimated mean of Y.
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
%   in_param.timebudget --- the time budget to do the two-stage estimation,
%   the default value is 100 seconds.
%
%   in_param.nbudget --- the sample budget to do the two-stage estimation,
%   the default value is 1e8.
%
%   in_param.npcmax --- number of elements in an array of optimal size to
%   calculate the mu, the default value is 1e6.
%
%   in_param.checked --- the status that the paramtered are checked.
%                        0   not checked
%                        1   checked by cubMC
%                        2   checked by meanMC
%
%   out_param_time_n_sigma_predict --- the estimated time to get n_sigma
%   samples of the random variable.
%
%   out_param_n_left_predict --- using the time left to predict the number
%   of samples left.
%
%   out_param.nmax --- the maximum sample budget to estimate mu, it comes
%   from both the sample budget and the time budget.
%
%   out_param.var --- the sample variance.
%
%   out_param.exit --- the state of program when exiting.
%                      0   Success.
%                      1   No enough samples to estimate the mean.                                
%                      2   Initial try out time costs more than 10% of time budget.                                 
%                      3   The estimated time for estimating variance is bigger
%                          than half of the time budget.
%                                      
%   out_param.kurtmax --- the upper bound on modified kurtosis.
%
%   out_param.n_mu --- the sample size that needed to estimate the mu.
%
%   out_param.n --- the total sample size needed to do the two stage
%   algorithm.
%
%   out_param.time --- the time eclipsed.
%
%  Guarantee
%
%  If the modified kurtosis of the random variable, Y, is less than the kmax,
%  which is defined in terms of the uncertainty, alpha, the sample size to
%  estimate variance, n_sigma, and the standard deviation inflation factor,
%  fudge, then the inequality
%
%  Pr(|mu-\hat{mu}| <= tol) >= 1-alpha 
%  
%  holds. Here mu is the true mean of Y, and \hat{mu} is the output
%  of MEANMC_G.
%
%  The cost of the two-stage algorithm also satisfies the inequality
%
%  Pr (N_tot <= N_up) >= 1-beta
%  
%  where N_tot is the total cost of samples, N_up is the upper bound on the
%  cost, which is roughly propotional to sigma^2/epsilon^2, beta is the
%  level of uncertainty on the cost. For details, please refer to [1].
%
%   Examples
%
%   Example 1: 
%   Calculate the mean of x^2 when x is uniformly distributed in
%   [0,1], with the absolute error tolerance = 1e-2.
%
%   >> in_param.abstol=1e-2; in_param.alpha = 0.01; Yrand=@(n) rand(n,1).^2;
%   >> mu=meanMC_g(Yrand,in_param) 
%   mu = 0.3***
%
%
%   Example 2: 
%   Using the same function as example 1, with the absolute error tolerance
%   1e-2.
%
%   >> mu=meanMC_g(Yrand,1e-2) 
%   mu = 0.3***
%
%
%   Example 3: 
%   Using the sample function as example 1, with the absolute error
%   tolerance 1e-2 and uncertainty 0.01.
%
%   >> mu=meanMC_g(Yrand,'abstol',1e-2,'alpha',0.01) 
%   mu = 0.3***
%
%
%   See also FUNAPPX_G, INTEGRAL_G, CUBMC_G
%
%   Reference
%   [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G.
%   W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to
%   appear, arXiv:1208.4318 [math.ST]

tstart = tic; %start the clock
[Yrand, in_param,out_param] = meanMC_g_param(varargin{:});

%out_param = in_param;%let the out_param contains all the in_param
n1 = 2;
Yrand(n1); %let it run once to load all the data. warm up the machine.
nsofar = n1;
ntry = 4; %initial try out sample size to get the time.
nsofar = nsofar + ntry;%the sample size has been used so far
tic;
Yrand(ntry);
timetry = toc; %count the time eclipsed to do ntry samples
while true
    timeleft = in_param.timebudget - toc(tstart);%update the time budget.
    out_param.n_left_predict = floor(timeleft*ntry/timetry);
    %get the estimated sample budget from time left
    out_param.time_n_sigma_predict = timetry/ntry * in_param.n_sigma;
    %get the estimated time using n_sigma samples.
    out_param.nmax = min(in_param.nbudget-nsofar,out_param.n_left_predict);
    %the max sample budget we could afford to do the following calculation
    if timetry > in_param.timebudget/10;
        %after intial try, found it has already used 10% of the time
        %budget, stop try.
        out_param.exit = 2; % exit the loop
        out_param = meanMC_g_err(out_param,tstart);% print the error message
        out_param.n_mu = out_param.nmax;
        break;
    elseif ntry >= ceil(in_param.n_sigma/100);
        %try out sample size is bigger than 1% of n_sigma that is used to
        %estimate the variance.
        if  out_param.time_n_sigma_predict > in_param.timebudget/2;
            % the estimated time using n_sigma samples is bigger than half
            % of the time budget, could not afford computing variance.
            % using all the sample left to compute the mean.
            out_param.exit = 3; % exit the loop
            out_param = meanMC_g_err(out_param,tstart);% print error message
            out_param.n_mu = out_param.nmax;
            break;
        else
            nsofar = nsofar+in_param.n_sigma;% the samples that have been used
            tic
            Yval = Yrand(in_param.n_sigma); % get the function values
            timetry_n_sigma = toc;
            % get the time for calculating n_sigma function values.
            out_param.var = var(Yval);% calculate the sample variance--stage 1
            sig0 = sqrt(out_param.var);% standard deviation
            sig0up = in_param.fudge.*sig0;% upper bound of the standard deviation
            alpha1 = 1-sqrt(1-in_param.alpha);% one side of the uncertainty
            out_param.kurtmax = (in_param.n_sigma-3)/(in_param.n_sigma-1) ...
                + ((alpha1*in_param.n_sigma)/(1-alpha1))*(1-1/in_param.fudge^2)^2;
            % get the upper bound on the modified kurtosis
            if sig0up == 0; % if the variance is zero, just take n_sigma samples
                out_param.n_mu = in_param.n_sigma;
                break;
            else
                toloversig = in_param.abstol/sig0up;
                % absolute error tolerance over sigma
                ncheb = ceil(1/(alpha1*toloversig.^2));
                % use Chebyshev inequality to estimate n
                A=18.1139;
                A1=0.3328;
                A2=0.429; % three constants in Berry-Esseen inequality
                M3upper=out_param.kurtmax^(3/4);%using Jensen inequality to
                % bound the third moment
                BEfun=@(logsqrtn)stdnormcdf(-exp(logsqrtn).*toloversig)...
                    +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
                    A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))- alpha1/2;
                % Berry-Esseen Inequality
                if BEfun(log(sqrt(in_param.n_sigma))) <= 0 || ...
                        ncheb <= in_param.n_sigma;
                    out_param.n_mu=in_param.n_sigma;
                    break;
                    %the Chebyshev n or the BE n is too small, just use
                    %param.n_sigma;
                else
                    logsqrtnCLT=log(stdnorminv(1-alpha1/2)/toloversig);
                    % get log of sqrt of n
                    out_param.n_mu=min(ncheb,ceil(exp(2*fzero(BEfun,logsqrtnCLT))));
                    % get the min n (used to estimate mu) by using cheb and BEfun
                    timeleft = in_param.timebudget-toc(tstart);% update the time left
                    out_param.n_left_predict = floor(timeleft*in_param.n_sigma...
                        /timetry_n_sigma);%using the time left to get the n left
                    out_param.nmax = min(in_param.nbudget-nsofar,...
                        out_param.n_left_predict);
                    % update the max n which will be used for estimating mu
                    if out_param.n_mu > out_param.nmax;
                        % if the sample size got from Chebyshev and BE fun is
                        % larger than nmax, print warning message and use nmax
                        out_param.exit=1; % exit the loop
                        meanMC_g_err(out_param,tstart); % print warning message
                        out_param.n_mu = out_param.nmax;% update n_mu
                        break;
                    else
                        break;
                    end
                end
            end
        end
    else
        multiplier = 5;
        ntry = multiplier*ntry; % boost the try out sample size five times
        nsofar=nsofar+ntry; % update the samples that have been used
        tic;
        Yrand(ntry);
        timetry = toc;% get the time for calclating ntry function values
    end
end
%%  Split The Param.n into columns
nopt=min(in_param.npcmax,out_param.n_mu);
% numbers of samples per loop step
nn=floor(out_param.n_mu/nopt); % number of loop steps
nremain=out_param.n_mu-nn*nopt;
% number of samples in last loop step
nloop=repmat(nopt,1,nn);
%vector of numbers of samples per loop step
if nremain>0; nloop=[nloop nremain]; nn=nn+1; end
sumY=0;
for iloop=1:nn %loops to save memory
    sumY=sumY+sum(Yrand(nloop(iloop)));
end
%%  Estimate mu
out_param.mu=sumY/out_param.n_mu; %calculate the mean
out_param.n=out_param.n_mu+nsofar;
%total number of samples used
mu=out_param.mu; %assign answer
out_param.time=toc(tstart); %elapsed time
end

function p = stdnormcdf(z)
% this function is to define cumulative distribution function (CDF) of the
% standard normal distribution.
p = 0.5*erfc(-z./sqrt(2));
% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
end

function x = stdnorminv(p)
% this function is the inverse function of CDF of standard normal distribution
x = -sqrt(2).*erfcinv(2*p);
end

function  [Yrand, in_param, out_param] = meanMC_g_param(varargin)

default.timebudget = 100;% default time budget
default.nbudget = 1e8; % default sample budget
default.abstol  = 1e-2;% default absolute error tolerance
default.n_sigma = 1e3;% default initial sample size n_sigma
default.fudge = 1.1;% default fudge factor
default.alpha = 0.01;% default uncertainty
default.npcmax = 1e6;% default n piece maximum
default.checked = 0;% default value for the paramters to be checked
if isempty(varargin)
    help meanMC_g
    warning('MATLAB:meanMC_g:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand = rand(n,1).^2.')
    Yrand = @(n) rand(n,1).^2;
    %give the error message
else
    Yrand = varargin{1};
end

validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if ~validvarargin
%if only have one input which is Yrand, use all the default parameters   
    in_param.abstol = default.abstol;
    in_param.alpha = default.alpha;
    in_param.n_sigma = default.n_sigma;
    in_param.fudge = default.fudge;
    in_param.timebudget = default.timebudget;
    in_param.nbudget = default.nbudget;
    in_param.npcmax = default.npcmax;
    in_param.checked = default.checked;
else
    p = inputParser;
    addRequired(p,'Yrand',@isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'n_sigma',default.n_sigma,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'timebudget',default.timebudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
        addOptional(p,'npcmax',default.npcmax,@isnumeric);
        addOptional(p,'checked',default.checked,@isnumeric);

    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'n_sigma',default.n_sigma,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'timebudget',default.timebudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
        addParamValue(p,'npcmax',default.npcmax,@isnumeric);
        addParamValue(p,'checked',default.checked,@isnumeric);
    end
    parse(p,Yrand,varargin{2:end})
    in_param = p.Results;
end

if in_param.checked==0
    if (in_param.abstol <= 0)% absolute error tolerance 
        warning('MATLAB:meanMC_g:abstolneg',...
            ['Absolute error tolerance should be greater than 0, ' ...
            'use the absolute value of the error tolerance'])
        in_param.abstol = abs(in_param.abstol);
    end
    if (in_param.alpha <= 0 ||in_param.alpha >= 1) % uncertainty
        warning('MATLAB:meanMC_g:alphanotin01',...
            ['the uncertainty should be less than 1 and bigger than 0, '...
            'use the default value.'])
        in_param.alpha = default.alpha;
    end
    if (~isposint(in_param.n_sigma)) % initial sample size should be an integer
        warning('MATLAB:meanMC_g:nsignotposint',...
            ['the number n_sigma should a positive integer, '...
            'take the absolute value and ceil.'])
        in_param.n_sigma = ceil(abs(in_param.n_sigma));
    end
    if (in_param.fudge <= 1) % standard deviation inflation factor/fudge factor
        warning('MATLAB:meanMC_g:fudgelessthan1',...
            ['the fudge factor should be bigger than 1, '...
            'use the default value.'])
        in_param.fudge = default.fudge;
    end
    if (in_param.timebudget <= 0) % time budget
        warning('MATLAB:meanMC_g:timebudgetlneg',...
            ['Time budget should be bigger than 0, '...
            'use the absolute value of time budget'])
        in_param.timebudget = abs(in_param.timebudget);
    end
    if (~isposint(in_param.nbudget)) % sample budget should be an integer
        warning('MATLAB:meanMC_g:nbudgetnotposint',...
            ['the number of sample budget should be a positive integer,'...
            'take the absolute value and ceil.'])
        in_param.nbudget = ceil(abs(in_param.nbudget));
    end
    if (~isposint(in_param.npcmax))
        % maxinum number of scalar values of x per vector should be a integer
        warning('MATLAB:meanMC_g:npcmaxnotposint',...
            ['the number of each piece of the samples should be' ...
            'a positive integer, take the absolute value and ceil.'])
        in_param.npcmax = ceil(abs(in_param.npcmax));
    end
in_param.checked = 2;
end
out_param = in_param;%let the out_param contains all the in_param
end

function [out_param,mu]=meanMC_g_err(out_param,tstart)
% Handles errors in meanMC_g and meanMC_g_param to give an exit with
%  information. 
%            out_param.exit = 0   success
%                             1   too many samples required
%                             2   too much try out time 
%                             3   too much time required to estimate
%                                 variance
if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
    case 1 % not enough samples to estimate the mean.
        warning('MATLAB:meanMC_g:maxreached',...
            ['tried to evalute at ' int2str(out_param.n_mu) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ' samples. Just use the maximum sample budget']);
        return
    case 2 % initial try out time costs more than 10% of time budget.
        warning('MATLAB:meanMC_g:initialtryoutbudgetreached',...
            ['initial try costs more than 10 percent '...
            'of time budget, stop try and return an answer '...
            'without guarantee.']);
        return
    case 3
        % the estimated time for estimating variance is bigger than half of
        % time budget.
        warning('MATLAB:meanMC_g:timebudgetreached',...
            ['the estimated time using n_sigma samples '...
            'is bigger than half of the time budget, '...
            'could not afford estimating variance, '...
            'use all the time left to estimate the mean.']);
        return
end
out_param.mu=NaN;
mu=out_param.mu;
if nargin>1; out_param.time=toc(tstart); end
end

