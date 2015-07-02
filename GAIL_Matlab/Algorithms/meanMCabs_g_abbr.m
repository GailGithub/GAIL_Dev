function [tmu,out_param]=meanMCabs_g_abbr(Yrand,abstol,alpha,nSig,fudge)
%MEANMCabs_g_abbr Monte Carlo method to estimate the mean of a random variable
%
%   tmu = MEANMCabs_g_abbr(Yrand,abstol,alpha,nSig,fudge) estimates the mean, mu, of a random variable Y to
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
%     percentage. The default value is 1%.
%
%     nSig --- the number of samples used to compute the sample variance
%
%     fudge --- the standard deviation invlation factor
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
if nargin < 5
    fudge = 1.2; %variance inflation factor
    if nargin < 4;
        nSig = 1e2; %number of samples to estimate variance
        if nargin < 3
            alpha = 0.01; %uncertainty
            if nargin < 2
                abstol = 0.01; %absolute error tolerance
                if nargin < 1
                    Yrand = @(n) rand(n,1).^2; %random number generator
                end
            end
        end
    end
end
nMax=1e8; %maximum number of samples allowed.
out_param.alpha = alpha; %save the input parameters to a structure
out_param.fudge = fudge;
out_param.nSig = nSig;
out_param.abstol = abstol;
tstart = tic; %start the clock
Yval = Yrand(nSig);% get samples to estimate variance
out_param.var = var(Yval); %calculate the sample variance--stage 1
sig0 = sqrt(out_param.var); %standard deviation
sig0up = out_param.fudge.*sig0; %upper bound on the standard deviation
alpha1 = 1-sqrt(1-out_param.alpha);% one side of the uncertainty
out_param.kurtmax = (out_param.nSig-3)/(out_param.nSig-1) ...
    + ((alpha1*out_param.nSig)/(1-alpha1))*(1-1/out_param.fudge^2)^2;
toloversig = out_param.abstol/sig0up;
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
logsqrtnCLT=log(stdnorminv(1-alpha1/2)/toloversig);
nmu=min(ncheb,ceil(exp(2*fzero(BEfun,logsqrtnCLT))));
% get the min n (used to estimate mu) by using cheb and BEfun
assert(nmu<nMax,['nmu = ' int2str(nmu) ', which is too big']) 
   %don't exceed sample budget
tmu=mean(Yrand(nmu)); %estimated mean
out_param.ntot=nSig+nmu; %total samples required
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


