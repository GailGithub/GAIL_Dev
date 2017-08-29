function [mu,out_param]=meanMCabs_g(varargin)
%MEANMCABS_G Monte Carlo method to estimate the mean of a random variable to
%within a specified absolute error tolerance with guaranteed confidence
%level 1-alpha.
%
%   mu = MEANMCABS_G(Yrand) estimates the mean of a random variable Y to within
%   a specified absolute error tolerance 1e-2 with guaranteed confidence
%   level 99%. Input Yrand is a function handle that accepts a positive
%   integer input n and returns an n x 1 vector of IID instances of the
%   random variable Y.
%
%   mu = MEANMCABS_G(Yrand,abstol,alpha,n_sigma,fudge,tbudget,nbudget,npcmax,checked)
%   estimates the mean of a random variable Y to within an specified absolute
%   error tolerance abstol with guaranteed confidence level 1-alpha. using
%   all ordered parsing inputs abstol, n_sigma, fudge, tbudget, nbudget,
%   npcmax and checked.
%
%   mu = MEANMCABS_G(Yrand,'abstol',abstol,'alpha',alpha,'n_sigma',n_sigma,...
%   'fudge',fudge,'tbudget',tbudget,'nbudget',nbudget,'npcmax',npcmax,...
%   'checked',checked) estimates the mean of a random variable Y to within a
%   specified absolute error tolerance abstol with guaranteed confidence
%   level 1-alpha. All the field-value pairs are optional and can be
%   supplied in different order.
%
%   mu = MEANMCABS_G(Yrand,in_param) estimates the mean of a random variable Y
%   to within a specified absolute error tolerance in_param.abstol with
%   guaranteed uncertainty within in_param.alpha. If a field is not
%   specified, the default value is used.
%
%   [mu, out_param] = MEANMCABS_G(Yrand,in_param) estimates the mean of a
%   random variable Y to within a specified absolute error tolerance with
%   the given parameters in_param and produce output parameters out_param.
%
%   Input Arguments
%
%     Yrand --- the function for generating IID instances of a random
%     variable Y whose mean we want to estimate. Y is often defined as a
%     function of some random variable X with a simple distribution.  For
%     example, if Y = X.^2 where X is a standard uniform random variable,
%     then one may define Yrand = @(n) rand(n,1).^2.
%
%     in_param.abstol --- the absolute error tolerance, default value is 1e-2.
% 
%     in_param.alpha --- the uncertainty, default value is 1%.
% 
%     in_param.n_sigma --- initial sample size for estimating the sample
%     variance, the default value is 1e4.
% 
%     in_param.fudge --- the standard deviation inflation factor, the default
%     value is 1.2.
% 
%     in_param.tbudget --- the time budget to do the two-stage estimation,
%     the default value is 100 seconds.
% 
%     in_param.nbudget --- the sample budget to do the two-stage estimation,
%     the default value is 1e8.
% 
%     in_param.npcmax --- number of elements in an array of optimal size to
%     calculate the mu, the default value is 1e6.
% 
%     in_param.checked --- the value corresponding to parameter checking status.
%                         0   not checked
%                         1   checked by meanMCabs_g
%
%   Output Arguments
%
%     mu --- the estimated mean of Y.
% 
%     out_param.time_n_sigma_predict --- the estimated time to get n_sigma
%     samples of the random variable.
%
%     out_param.n_left_predict --- using the time left to predict the number
%     of samples left.
% 
%     out_param.nmax --- the maximum sample budget to estimate mu, it comes
%     from both the sample budget and the time budget.
% 
%     out_param.var --- the sample variance.
%
%     out_param.exit --- the state of program when exiting.
%                       0   Success.
%                       1   No enough samples to estimate the mean.                                
%                       2   Initial try out time costs more than 10% of time budget.                                 
%                       3   The estimated time for estimating variance is bigger
%                           than half of the time budget.
% 
%     out_param.kurtmax --- the upper bound on modified kurtosis.
% 
%     out_param.n_mu --- the sample size that needed to estimate the mu.
% 
%     out_param.n --- the total sample size needed to do the two stage estimation.
% 
%     out_param.time --- the time elapsed.
% 
%  Guarantee
% 
% If the modified kurtosis of the random variable, Y, is less than the kurtmax,
% which is defined in terms of the uncertainty, alpha, the sample size to
% estimate variance, n_sigma, and the standard deviation inflation factor,
% fudge, then the inequality
% 
% Pr(|mu-\hat{mu}| <= abstol) >= 1-alpha 
% 
% holds. Here mu is the true mean of Y, and \hat{mu} is the output
% of MEANMCABS_G.
% 
% The cost of the two-stage algorithm also satisfies the inequality
% 
% Pr (N_tot <= N_up) >= 1-beta
% 
% where N_tot is the total cost of samples, N_up is the upper bound on the
% cost, which is roughly proportional to sigma^2/abstol^2, beta is the
% level of uncertainty on the cost. For details, please refer to [1].
% 
%  Examples
% 
% Example 1: 
% Calculate the mean of x^2 when x is uniformly distributed in
% [0 1], with the absolute error tolerance = 1e-2.
% 
% >> in_param.abstol=1e-2; in_param.alpha = 0.01; Yrand=@(n) rand(n,1).^2;
% >> mu=meanMCabs_g(Yrand,in_param) 
% mu = 0.3***
% 
% 
% Example 2: 
% Using the same function as example 1, with the absolute error tolerance
% 1e-2.
% 
% >> mu=meanMCabs_g(Yrand,1e-2) 
% mu = 0.3***
% 
% 
% Example 3: 
% Using the sample function as example 1, with the absolute error
% tolerance 1e-2 and uncertainty 0.01.
% 
% >> mu=meanMCabs_g(Yrand,'abstol',1e-2,'alpha',0.01) 
% mu = 0.3***
% 
% 
%   See also FUNAPPX_G, INTEGRAL_G, CUBMCABS_G
%
%  References
%
%   [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
%   Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
%   arXiv:1208.4318 [math.ST]
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
%   Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
%   1.3.0)" [MATLAB Software], 2014. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above paper and software.

tstart = tic; %start the clock
[Yrand, out_param] = meanMCabs_g_param(varargin{:});

n1 = 2;
Yrand(n1); %let it run once to load all the data. warm up the machine.
nsofar = n1;
ntry = 4; %initial try out sample size to get the time.
nsofar = nsofar + ntry;%the sample size has been used so far
tic;
Yrand(ntry);
timetry = toc; %count the time elapsed to do ntry samples
out_param.exit = 0;%assume success until finding error
while true
    timeleft = out_param.tbudget - toc(tstart);%update the time budget.
    out_param.n_left_predict = floor(timeleft*ntry/timetry);
    %get the estimated sample budget from time left
    out_param.time_n_sigma_predict = timetry/ntry * out_param.n_sigma;
    %get the estimated time using n_sigma samples.
    out_param.nmax = max(min(out_param.nbudget-nsofar,out_param.n_left_predict),1);
    %the max sample budget we could afford to do the following calculation
    if timetry > out_param.tbudget/10;
        %after intial try, found it has already used 10% of the time
        %budget, stop try.
        out_param.exit = 2; % exit the loop
        out_param = meanMCabs_g_err(out_param);% print the error message
        out_param.n_mu = out_param.nmax;
        break;
    elseif ntry >= ceil(out_param.n_sigma/100);
        %try out sample size is bigger than 1% of n_sigma that is used to
        %estimate the variance.
        if  out_param.time_n_sigma_predict > out_param.tbudget/2;
            % the estimated time using n_sigma samples is bigger than half
            % of the time budget, could not afford computing variance.
            % using all the sample left to compute the mean.
            out_param.exit = 3; % exit the loop
            out_param = meanMCabs_g_err(out_param);% print error message
            out_param.n_mu = out_param.nmax;
            break;
        else
            nsofar = nsofar+out_param.n_sigma;% the samples that have been used
            tic
            Yval = Yrand(out_param.n_sigma); % get the function values
            timetry_n_sigma = toc;
            % get the time for calculating n_sigma function values.
            out_param.var = var(Yval);% calculate the sample variance--stage 1
            sig0 = sqrt(out_param.var);% standard deviation
            sig0up = out_param.fudge.*sig0;% upper bound of the standard deviation
            alpha1 = 1-sqrt(1-out_param.alpha);% one side of the uncertainty
            out_param.kurtmax = (out_param.n_sigma-3)/(out_param.n_sigma-1) ...
                + ((alpha1*out_param.n_sigma)/(1-alpha1))*(1-1/out_param.fudge^2)^2;
            % get the upper bound on the modified kurtosis
            if sig0up == 0; %#ok<BDSCI> % if the variance is zero, just take n_sigma samples
                out_param.n_mu = out_param.n_sigma;
                break;
            else
                toloversig = out_param.abstol/sig0up;
                % absolute error tolerance over sigma
                ncheb = ceil(1/(alpha1*toloversig.^2));
                % use Chebyshev inequality to estimate n
%                 A=18.1139;
%                 A1=0.3328;
%                 A2=0.429; % three constants in Berry-Esseen inequality
%                 M3upper=out_param.kurtmax^(3/4);%using Jensen inequality to
%                 % bound the third moment
%                 BEfun=@(logsqrtn)stdnormcdf(-exp(logsqrtn).*toloversig)...
%                     +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
%                     A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))- alpha1/2;
%                 % Berry-Esseen Inequality
                A=18.1139;
                A1=0.3322;
                A2=0.429;
                A3= 0.3031;
                A4= 0.646;
                A5= 0.469; % Six constants in Berry-Esseen inequality
                M3upper=out_param.kurtmax^(3/4);%using Jensen inequality to
                % bound the third moment
                BEfun=@(logsqrtn)stdnormcdf(-exp(logsqrtn).*toloversig)...
                    +exp(-logsqrtn).*min([A1*(M3upper+A2),A3*(M3upper+A4),A5*M3upper, ...
                    A*M3upper./(1+(exp(logsqrtn).*toloversig).^3)])- alpha1/2;
                % Berry-Esseen Inequality               
                if BEfun(log(sqrt(out_param.n_sigma))) <= 0 || ...
                        ncheb <= out_param.n_sigma;
                    out_param.n_mu=out_param.n_sigma;
                    break;
                    %the Chebyshev n or the BE n is too small, just use
                    %param.n_sigma;
                else
                    logsqrtnCLT=log(stdnorminv(1-alpha1/2)/toloversig);
                    % get log of sqrt of n
                    out_param.n_mu=min(ncheb,ceil(exp(2*fzero(BEfun,logsqrtnCLT))));
                    % get the min n (used to estimate mu) by using cheb and BEfun
                    timeleft = out_param.tbudget-toc(tstart);% update the time left
                    out_param.n_left_predict = floor(timeleft*out_param.n_sigma...
                        /timetry_n_sigma);%using the time left to get the n left
                    out_param.nmax = max(min(out_param.nbudget-nsofar,...
                        out_param.n_left_predict),1);
                    % update the max n which will be used for estimating mu
                    if out_param.n_mu > out_param.nmax;
                        % if the sample size got from Chebyshev and BE fun is
                        % larger than nmax, print warning message and use nmax
                        out_param.exit=1; % exit the loop
                        meanMCabs_g_err(out_param); % print warning message
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
        ntry = multiplier*ntry; % boost the try out sample size at multiplier times
        nsofar=nsofar+ntry; % update the samples that have been used
        tic;
        Yrand(ntry);
        timetry = toc;% get the time for calculating ntry function values
    end
end
%%  Split The Param.n into columns
nopt=min(out_param.npcmax,out_param.n_mu);
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

function  [Yrand,out_param] = meanMCabs_g_param(varargin)
default.tbudget = 100;% default time budget
default.nbudget = 1e8; % default sample budget
default.abstol  = 1e-2;% default absolute error tolerance
default.n_sigma = 1e4;% default initial sample size n_sigma
default.fudge = 1.2;% default fudge factor
default.alpha = 0.01;% default uncertainty
default.npcmax = 1e6;% default n piece maximum
default.checked = 0;% default value of parameter checking status

if isempty(varargin)
    help meanMCabs_g
    warning('GAIL:meanMCabs_g:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand = rand(n,1).^2.')
    Yrand = @(n) rand(n,1).^2;
    %give the error message
else
    Yrand = varargin{1};
    if max(size(Yrand(5)))~=5 || min(size(Yrand(5)))~=1
        % if the input is not a length n Vector, print warning message
        warning('GAIL:meanMCabs_g:yrandnotlengthN',...
            ['Yrand should be a random variable vector of length n, '...
            'but not an integrand or a matrix'])
    end
end

validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if ~validvarargin
%if only have one input which is Yrand, use all the default parameters   
    out_param.abstol = default.abstol;
    out_param.alpha = default.alpha;
    out_param.n_sigma = default.n_sigma;
    out_param.fudge = default.fudge;
    out_param.tbudget = default.tbudget;
    out_param.nbudget = default.nbudget;
    out_param.npcmax = default.npcmax;
    out_param.checked = default.checked;
else
    p = inputParser;
    addRequired(p,'Yrand',@gail.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'n_sigma',default.n_sigma,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
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
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
        addParamValue(p,'npcmax',default.npcmax,@isnumeric);
        addParamValue(p,'checked',default.checked,@isnumeric);
    end
    parse(p,Yrand,varargin{2:end})
    out_param = p.Results;
end

if out_param.checked==0
    if (out_param.abstol <= 0)% absolute error tolerance 
        warning('GAIL:meanMCabs_g:abstolneg',...
            ['Absolute error tolerance should be greater than 0, ' ...
            'use the absolute value of the error tolerance'])
        out_param.abstol = abs(out_param.abstol);
    end
    if (out_param.alpha <= 0 ||out_param.alpha >= 1) % uncertainty
        warning('GAIL:meanMCabs_g:alphanotin01',...
            ['the uncertainty should be less than 1 and bigger than 0, '...
            'use the default value.'])
        out_param.alpha = default.alpha;
    end
    if (~gail.isposint(out_param.n_sigma)) % initial sample size should be an integer
        warning('GAIL:meanMCabs_g:nsignotposint',...
            ['the number n_sigma should a positive integer, '...
            'take the absolute value and ceil.'])
        out_param.n_sigma = ceil(abs(out_param.n_sigma));
    end
    if (out_param.fudge <= 1) % standard deviation inflation factor/fudge factor
        warning('GAIL:meanMCabs_g:fudgelessthan1',...
            ['the fudge factor should be bigger than 1, '...
            'use the default value.'])
        out_param.fudge = default.fudge;
    end
    if (out_param.tbudget <= 0) % time budget
        warning('GAIL:meanMCabs_g:timebudgetlneg',...
            ['Time budget should be bigger than 0, '...
            'use the absolute value of time budget'])
        out_param.tbudget = abs(out_param.tbudget);
    end
    if (~gail.isposint(out_param.nbudget)) % sample budget should be an integer
        warning('GAIL:meanMCabs_g:nbudgetnotposint',...
            ['the number of sample budget should be a positive integer,'...
            'take the absolute value and ceil.'])
        out_param.nbudget = ceil(abs(out_param.nbudget));
    end
    if (~gail.isposint(out_param.npcmax))
        % maximum number of scalar values of x per vector should be a integer
        warning('GAIL:meanMCabs_g:npcmaxnotposint',...
            ['the number of each piece of the samples should be' ...
            'a positive integer, take the absolute value and ceil.'])
        out_param.npcmax = ceil(abs(out_param.npcmax));
    end
out_param.checked = 1;
end
end

function out_param = meanMCabs_g_err(out_param)
% Handles errors in meanMCabs_g and meanMCabs_g_param to give an exit with
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
        warning('GAIL:meanMCabs_g:maxreached',...
            ['tried to evaluate at ' int2str(out_param.n_mu) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ' samples. Just use the maximum sample budget.']);
        return
    case 2 % initial try out time costs more than 10% of time budget.
        warning('GAIL:meanMCabs_g:initialtryoutbudgetreached',...
            ['initial try costs more than 10 percent '...
            'of time budget, stop try and return an answer '...
            'without guarantee.']);
        return
    case 3
        % the estimated time for estimating variance is bigger than half of
        % time budget.
        warning('GAIL:meanMCabs_g:timebudgetreached',...
            ['the estimated time using n_sigma samples '...
            'is bigger than half of the time budget, '...
            'could not afford estimating variance, '...
            'use all the time left to estimate the mean.']);
        return
end
end

