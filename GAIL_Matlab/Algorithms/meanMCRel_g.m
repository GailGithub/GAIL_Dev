function [mu,out_param]=meanMCRel_g(varargin)
% MEANMCRel_G Monte Carlo method to estimate the mean of a random variable to
% within a specified generalized error tolerance with guaranteed confidence
% level 1-alpha.
%
% mu = MEANMCRel_G(Yrand) estimates the mean of a random variable Y to within
% a specified generalized error tolerance 1e-2 with guaranteed confidence
% level 99%. Input Yrand is a function handle that accepts a positive
% integer input n and returns an n x 1 vector of IID instances of the
% random variable Y.
%
% mu =
% MEANMCRel_G(Yrand,abstol,alpha,n_sigma,fudge,tbudget,nbudget,npcmax,checked)
% estimates the mean of a random variable Y to within an specified generalized
% error tolerance abstol with guaranteed confidence level 1-alpha. using
% all ordered parsing inputs abstol, n_sigma, fudge, tbudget, nbudget,
% npcmax and checked.
%
% mu = MEANMCRel_G(Yrand,'abstol',abstol,'alpha',alpha,'n_sigma',n_sigma,...
% 'fudge',fudge,'tbudget',tbudget,'nbudget',nbudget,'npcmax',npcmax,...
% 'checked',checked) estimates the mean of a random variable Y to within a
% specified generalized error tolerance abstol with guaranteed confidence
% level 1-alpha. All the field-value pairs are optional and can be
% supplied in different order.
%
% mu = MEANMCRel_G(Yrand,in_param) estimates the mean of a random variable Y
% to within a specified generalized error tolerance in_param.abstol with
% guaranteed uncertainty within in_param.alpha. If a field is not
% specified, the default value is used.
%
% [mu, out_param] = MEANMCRel_G(Yrand,in_param) estimates the mean of a
% random variable Y to within a specified generalized error tolerance with
% the given parameters in_param and produce output parameters out_param.
%
% Yrand --- the function for generating IID instances of a random
% variable Y whose mean we want to estimate. Y is often defined as a
% function of some random variable X with a simple distribution.  For
% example, if Y = X.^2 where X is a standard uniform random variable,
% then one may define Yrand = @(n) rand(n,1).^2.
%
% mu --- the estimated mean of Y.
%
% in_param.abstol --- the absolute error tolerance, default value is 1e-2.
%
% in_param.reltol --- the relative error tolerance, default value is 1e-1.
%
% in_param.alpha --- the uncertainty, default value is 1%.
%
% in_param.nSig --- initial sample size for estimating the sample
% variance, the default value is 1e3.
%
% in_param.n1 --- initial sample size for estimating the sample
% mean, the default value is 1e4.
%
% in_param.tbudget --- the time budget to do the two-stage estimation,
% the default value is 100 seconds.
%
% in_param.nbudget --- the sample budget to do the two-stage estimation,
% the default value is 1e8.
%
% in_param.checked --- the value corresponding to parameter checking status.
%                     0   not checked
%                     1   checked by cubMC_g
%                     2   checked by meanMC_g
%
% out_param.nmax --- the maximum sample budget to estimate mu, it comes
% from both the sample budget and the time budget.
%
% out_param.var --- the sample variance.
%
% out_param.exit --- the state of program when exiting.
%                  0   Success.
%                  1   No enough samples to estimate the mean.
%                  2   Initial try out time costs more than 10% of time budget.
%                  3   The estimated time for estimating variance is bigger
%                      than half of the time budget.
%
% out_param.kurtmax --- the upper bound on modified kurtosis.
%
% out_param.n_mu --- the sample size that needed to estimate the mu.
%
% out_param.n --- the total sample size used to do the two stage estimation.
%
% out_param.time --- the time elapsed.
%
% Guarantee
%
% If the modified kurtosis of the random variable, Y, is less than the kurtmax,
% which is defined in terms of the uncertainty, alpha, the sample size to
% estimate variance, n_sigma, and the standard deviation inflation factor,
% fudge, then the inequality
%
% Pr(|mu-\hat{mu}| <= abstol) >= 1-alpha
%
% holds. Here mu is the true mean of Y, and \hat{mu} is the output
% of MEANMC_G.
%
% The cost of the two-stage algorithm also satisfies the inequality
%
% Pr (N_tot <= N_up) >= 1-beta
%
% where N_tot is the total cost of samples, N_up is the upper bound on the
% cost, which is roughly proportional to sigma^2/abstol^2, beta is the
% level of uncertainty on the cost. For details, please refer to [1].
%
% Examples
%
% Example 1:
% Calculate the mean of x^2 when x is uniformly distributed in
% [0 1], with the absolute error tolerance = 1e-2.
%
% >> in_param.abstol=1e-2; in_param.alpha = 0.01; Yrand=@(n) rand(n,1).^2;
% >> mu=meanMCRel_g(Yrand,in_param)
% mu = 0.3***
%
%
% Example 2:
% Using the same function as example 1, with the absolute error tolerance
% 1e-2.
%
% >> mu=meanMCRel_g(Yrand,1e-2,1e-2)
% mu = 0.33***
%
%
% Example 3:
% Using the sample function as example 1, with the absolute error
% tolerance 1e-2 and uncertainty 0.01.
%
% >> mu=meanMCRel_g(Yrand,'abstol',1e-2,'alpha',0.01)
% mu = 0.3***
%
%
% See also FUNAPPX_G, INTEGRAL_G, CUBMC_G
%
% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
% conservative fixed width confidence intervals via Monte Carlo sampling,
% Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
% Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
% arXiv:1208.4318 [math.ST]
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
% Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
% 1.3.0)" [MATLAB Software], 2014. Available from
% http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.

tstart = tic; %start the clock
[Yrand, out_param] = meanMC_g_param(varargin{:});
npcmax=1e6;
fudge=1.1;
n1 = 2;
Yrand(n1); %let it run once to load all the data. warm up the machine.
nsofar = n1;
ntry = 4; %initial try out sample size to get the time.
nsofar = nsofar + ntry;%the sample size has been used so far
tic;
Yrand(ntry);
timetry = toc; %count the time elapsed to do ntry samples
while true
    out_param.nmax = getNleft(out_param.tbudget,tstart...
        ,ntry,timetry,nsofar,out_param.nbudget);
    %update the nmax after initinal try
    t_sig_predict = timetry/ntry * out_param.nSig;
    %get the estimated time using n_sigma samples.
    if timetry > out_param.tbudget/10;
        %after intial try, found it has already used 10% of the time
        %budget, stop try.
        out_param.exit = 2; % exit the loop
        out_param = meanMC_g_err(out_param);% print the error message
        out_param.n_mu = out_param.nmax;
        mu = SplitColumnMean(Yrand,out_param.n_mu,npcmax);
        %calculate the mean without guarantee using nmax
        out_param.ntot = out_param.n_mu+nsofar;%update the total sample used
        break;
    elseif ntry >= ceil(out_param.nSig/100);
        %try out sample size is bigger than 1% of nSig that is used to
        %estimate the variance.
        if  t_sig_predict > out_param.tbudget/2;
            % the estimated time using nSig samples is bigger than half
            % of the time budget, could not afford computing variance.
            % using all the sample left to compute the mean.
            out_param.exit = 3; % exit the loop
            out_param = meanMC_g_err(out_param);% print error message
            out_param.n_mu = out_param.nmax;
            mu = SplitColumnMean(Yrand,out_param.n_mu,npcmax);
            out_param.ntot = out_param.n_mu+nsofar;
            break;
        else
            tic;
            Yval = Yrand(out_param.nSig); % get the function values
            t_sig = toc;
            % get the time for calculating nSig function values.
            nsofar = nsofar+out_param.nSig;
            % update the samples that have been used
            out_param.var = var(Yval);% calculate the sample variance--stage 1
            sig0 = sqrt(out_param.var);% standard deviation
            sig0up = fudge.*sig0;% upper bound on standard deviation
            alpha_sig = out_param.alpha/2;% the uncertainty for variance estimation
            out_param.alphai(1) = (out_param.alpha-alpha_sig)/(2*(1-alpha_sig));
            out_param.kurtmax = (out_param.nSig-3)/(out_param.nSig-1) ...
                + ((alpha_sig*out_param.nSig)/(1-alpha_sig))...
                *(1-1/fudge^2)^2;
            % get the upper bound on the modified kurtosis
            ncbinv = N_CB_inv(out_param.n1,out_param.alphai(1),out_param.kurtmax);
            out_param.tol(1) = sig0up*ncbinv;
            % the width of initial confidence interval for the mean
            out_param.nmax = getNleft(out_param.tbudget,...
                tstart,out_param.nSig,t_sig,nsofar,out_param.nbudget);
            i=1;
            while true
                out_param.n(1) = out_param.n1;
                if out_param.n(i) > out_param.nmax;
                    % if the sample size used for initial estimation is
                    % larger than nmax, print warning message and use nmax
                    out_param.nexceed = out_param.n(i);
                    out_param.exit=1; % exit the loop
                    meanMC_g_err(out_param); % print warning message
                    out_param.n(i) = out_param.nmax;% update n1
                    mu = SplitColumnMean(Yrand,out_param.n(i),npcmax);
                    out_param.ntot = out_param.n(i)+nsofar;
                    break;
                end
                out_param.alphai(i+1) = (out_param.alpha-alpha_sig)/...
                    (1-alpha_sig)*2.^(-i-1);
                out_param.mu(i) = SplitColumnMean(Yrand,out_param.n(i),npcmax);
                nsofar = out_param.n(i)+nsofar;
                out_param.nmax = out_param.nmax-out_param.n(i);
                theta = 2;
                index = 1;
                out_param.DeltaPlus(i) = (tol_fun(out_param.abstol,out_param.reltol,...
                    theta,out_param.mu(i) - out_param.tol(i),index)...
                    +tol_fun(out_param.abstol,out_param.reltol,...
                    theta,out_param.mu(i) + out_param.tol(i),index))/2;
                if out_param.DeltaPlus(i) >= out_param.tol(i)
                    T = i;
                    mu = out_param.mu(T);
                    break;
                else
                    deltat=0.7;
                    deltah=0.5;
                    delta=0.3;
                    out_param.tol(i+1) = max(min(out_param.DeltaPlus(i)*deltat, ...
                        deltah*out_param.tol(i)),delta*out_param.tol(i));
                    toloversig = out_param.tol(i+1)/sig0up;
                    out_param.n(i+1) = N_CB(toloversig,out_param.alphai(i+1)...
                        ,out_param.kurtmax);
                    i=i+1;
                end
            end
            break;
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
out_param.ntot=sum(out_param.n)+out_param.nSig;
out_param.time=toc(tstart); %elapsed time
end

function mu =SplitColumnMean(RV,n,npcmax)
%%  Split The Param.n into columns and calculate the mean
nopt=min(npcmax,n);
% numbers of samples per loop step
nn=floor(n/nopt); % number of loop steps
nremain=n-nn*nopt;
% number of samples in last loop step
nloop=repmat(nopt,1,nn);
%vector of numbers of samples per loop step
if nremain>0; nloop=[nloop nremain]; nn=nn+1; end
sumY=0;
for iloop=1:nn %loops to save memory
    sumY=sumY+sum(RV(nloop(iloop)));
end
%%  Estimate mu
mu=sumY/n; %calculate the mean
end

function ncb = N_CB(toloversig,alpha,kmax)
ncheb = ceil(1/(toloversig^2*alpha));
A=18.1139;
A1=0.3328;
A2=0.429; % three constants in Berry-Esseen inequality
M3upper = kmax^(3/4);
BEfun2=@(logsqrtn)stdnormcdf(-exp(logsqrtn).*toloversig)...
    +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
    A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))- ...
    alpha/2;
logsqrtnCLT=log(stdnorminv(1-alpha/2)/toloversig);
% get log of sqrt of n
nbe=ceil(exp(2*fzero(BEfun2,logsqrtnCLT)));
ncb = min(ncheb,nbe);
end

function ncbinv = N_CB_inv(n1,alpha1,kmax)
NCheb_inv = 1/sqrt(n1*alpha1);
% use Chebyshev inequality
A=18.1139;
A1=0.3328;
A2=0.429; % three constants in Berry-Esseen inequality
M3upper=kmax^(3/4);%using Jensen inequality to
% bound the third moment
BEfun=@(logsqrtb)stdnormcdf(n1.*logsqrtb)...
    +min(A1*(M3upper+A2), ...
    A*M3upper./(1+(sqrt(n1).*logsqrtb).^3))/sqrt(n1)...
    - alpha1/2;
% Berry-Esseen Inequality
logsqrtb_clt=log(sqrt(stdnorminv(1-alpha1/2)/sqrt(n1)));
NBE_inv = exp(2*fzero(BEfun,logsqrtb_clt));
ncbinv = min(NCheb_inv,NBE_inv);
end

function eps = tol_fun(abstol,reltol,theta,mu,index)
switch index
    case 1 % the theta case
        %theta=0---absolute error
        %theta=1---relative error
        eps  = theta*abstol+(1-theta)*reltol;
    case 2 % the max case
        eps  = max(abstol,reltol*abs(mu));
end
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

function  [Yrand,out_param] = meanMC_g_param(varargin)

default.abstol  = 1e-2;% default absolute error tolerance
default.reltol = 1e-1;% default relative error tolerance
default.nSig = 1e2;% default initial sample size n_sigma for variance estimation
default.n1 = 1e3; % default initial sample size n1 for mean estimation
default.alpha = 0.01;% default uncertainty
default.tbudget = 200;% default time budget
default.nbudget = 1e9; % default sample budget
default.checked = 0;% default value of parameter checking status

if isempty(varargin)
    help meanMCRel_g
    warning('MATLAB:meanMCRel_g:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand = rand(n,1).^2.')
    Yrand = @(n) rand(n,1).^2;
    %give the error message
else
    Yrand = varargin{1};
    if max(size(Yrand(5)))~=5 || min(size(Yrand(5)))~=1
        % if the input is not a length n Vector, print warning message
        warning('MATLAB:meanMCRel_g:yrandnotlengthN',...
            ['Yrand should be a random variable vector of length n, '...
            'but not an integrand or a matrix'])
        Yrand = @(n) rand(n,1).^2;
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

    out_param.abstol = default.abstol;% default absolute error tolerance
    out_param.reltol = default.reltol; % default relative error tolerance
    out_param.n1 = default.n1;
    out_param.nSig = default.nSig;
    out_param.alpha = default.alpha;
    out_param.tbudget = default.tbudget;
    out_param.nbudget = default.nbudget;
    out_param.checked = default.checked;
else
    p = inputParser;
    addRequired(p,'Yrand',@GAIL_Internal.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.

        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'nSig',default.nSig,@isnumeric);
        addOptional(p,'n1',default.n1,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
        addOptional(p,'checked',default.checked,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end

        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'nSig',default.nSig,@isnumeric);
        addParamValue(p,'n1',default.n1,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);        
        addParamValue(p,'checked',default.checked,@isnumeric);
    end
    parse(p,Yrand,varargin{2:end})
    out_param = p.Results;
end

if out_param.checked==0
    if (out_param.abstol <= 0)
        % absolute error tolerance
        warning('MATLAB:meanMCRel_g:abstolneg',...
            ['Absolute error tolerance should be greater than 0, ' ...
            'use the absolute value of the error tolerance'])
        out_param.abstol = abs(out_param.abstol);
    end
     if (out_param.reltol <= 0 || out_param.reltol >= 1)
         % relative error tolerance
        warning('MATLAB:meanMCRel_g:reltolneg',...
            ['Relative error tolerance should be between 0 and 1, ' ...
            'use the default value of the error tolerance'])
        out_param.abstol = abs(out_param.abstol);
    end   
    if (out_param.alpha <= 0 ||out_param.alpha >= 1) % uncertainty
        warning('MATLAB:meanMCRel_g:alphanotin01',...
            ['the uncertainty should be between 0 and 1, '...
            'use the default value.'])
        out_param.alpha = default.alpha;
    end
    if (~GAIL_Internal.isposint(out_param.nSig)) % initial sample size should be an integer
        warning('MATLAB:meanMCRel_g:nsignotposint',...
            ['the number n_sigma should a positive integer, '...
            'take the absolute value and ceil.'])
        out_param.nSig = ceil(abs(out_param.nSig));
    end
    if (out_param.tbudget <= 0) % time budget should be positive
        warning('MATLAB:meanMCRel_g:timebudgetlneg',...
            ['Time budget should be bigger than 0, '...
            'use the absolute value of time budget'])
        out_param.tbudget = abs(out_param.tbudget);
    end
    if (~GAIL_Internal.isposint(out_param.nbudget)) % sample budget should be an integer
        warning('MATLAB:meanMCRel_g:nbudgetnotposint',...
            ['the number of sample budget should be a positive integer,'...
            'take the absolute value and ceil.'])
        out_param.nbudget = ceil(abs(out_param.nbudget));
    end
    out_param.checked = 2;
end
end

function out_param = meanMC_g_err(out_param)
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
        warning('MATLAB:meanMCRel_g:maxreached',...
            ['tried to evaluate at ' int2str(out_param.nexceed) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ' samples. Just use the maximum sample budget.']);
        return
    case 2 % initial try out time costs more than 10% of time budget.
        warning('MATLAB:meanMCRel_g:initialtryoutbudgetreached',...
            ['initial try costs more than 10 percent '...
            'of time budget, stop try and return an answer '...
            'without guarantee.']);
        return
    case 3
        % the estimated time for estimating variance is bigger than half of
        % time budget.
        warning('MATLAB:meanMCRel_g:timebudgetreached',...
            ['the estimated time using n_sigma samples '...
            'is bigger than half of the time budget, '...
            'could not afford estimating variance, '...
            'use all the time left to estimate the mean.']);
        return
end
end

function nmax= getNleft(tbudget,tstart,ntry,ttry,nsofar,nbudget)
            timeleft = tbudget-toc(tstart);% update the time left
            nleft = floor(timeleft*ntry/ttry);
            nmax = max(min(nbudget-nsofar,nleft),1);
end