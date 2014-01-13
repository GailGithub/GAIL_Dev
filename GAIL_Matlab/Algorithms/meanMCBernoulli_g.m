function [mu,out_param,in_param]=meanMCBernoulli_g(varargin)
% MEANMCBERNOULLI_G Monte Carlo method to estimate the mean of a Bernoullii
% random variable to within a specific absolute error tolerance with
% guaranteed uncertainty within alpha.
%
%   mu = MEANMCBERNOULLI_G(Yrand) estimates the mean of a Bernoulli random
%   variable Y to within a specified absolute error tolerance 1e-2 with
%   guaranteed uncertainty within 1%. Input Yrand is a function handle that
%   accepts a positive integer input n and returns an n x 1 vector of IID
%   instances of the Bernoulli random variable Y.
%
%   mu =
%   MEANMCBERNOULLI_G(Yrand,abstol,alpha,n_sigma,fudge,timebudget,nbudget,npcmax)
%   estimates the mean of a Bernoulli random variable Y to within an
%   absolute error tolerance abstol with guaranteed uncertainty within
%   alpha using all ordered parsing inputs n_sigma, fudge, timebudget,
%   nbudget and npcmax.
%
%   mu =
%   MEANMCBERNOULLI_G(Yrand,'abstol',abstol,'alpha',alpha,'n_sigma',n_sigma,
%   'fudge',fudge,'timebudget',timebudget,'nbudget',nbudget,'npcmax',npcmax)
%   estimates the mean of a Bernoulli random variable Y to within a
%   specified absolute error tolerance abstol with guaranteed untertainty
%   within alpha. All the field-value pairs are optional and can be
%   supplied in different order.
%
%   mu = MEANMCBERNOULLI_G(Yrand,in_param) estimates the mean of a
%   Bernoulli random variable Y to within a specified absolute error
%   tolerance in_param.abstol with guaranteed uncertainty within
%   in_param.alpha. If a field is not specified, the default value is used.
%
%   [mu, out_param] = MEANMCBERNOULLI_G(Yrand,in_param) estimates the mean
%   of a Bernoulli random variable Y to within a specified absolute error
%   tolerance with the given parameters in_param and output parameters
%   out_param.
%
%
%   Yrand --- the function for generating IID instances of a Bernoulli random
%   variable Y whose mean we want to estimate. 
%
%   mu --- the estimated mean of Y.
%
%   in_param.abstol --- the absolute error tolerance, default value is 1e-2.
%
%   in_param.alpha --- the uncertainty, default value is 1%.
%
%   in_param.npcmax --- number of elements in an array of optimal size to
%   calculate the mu, the default value is 1e6.
%
%   in_param.nmax --- the sample budget.
%
%   out_param.n --- the total sample size used.
%
%   out_param.time --- the time eclipsed.
%
%   Examples
%
%   Example 1: 
%   Calculate the mean of a bernoulli random variable with true mu=0.55,with
%   error tolerance 1e-3 and uncertainty 0.01.
%
%   >> in_param.abstol=1e-3; in_param.alpha = 0.01; p=1/90;Yrand=@(n) binornd(1,p,n,1);
%   >> mu=meanMCBernoulli_g(Yrand,in_param) 
%   mu = 0.01***
%
%
%   Example 2: 
%   Using the same function as example 1, with the absolute error tolerance
%   1e-2.
%
%   >> mu=meanMCBernoulli_g(Yrand,1e-4) 
%   mu = 0.011***
%
%
%   Example 3: 
%   Using the sample function as example 1, with the absolute error
%   tolerance 1e-4 and uncertainty 0.005.
%
%   >> mu=meanMCBernoulli_g(Yrand,'abstol',1e-4,'alpha',0.005) 
%   mu = 0.0111***
%
%
%   See also FUNAPPX_G, INTEGRAL_G, CUBMC_G, MEANMC_G
%
%   Reference
%   [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G.
%   W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to
%   appear, arXiv:1208.4318 [math.ST]
%
tstart = tic; %start the clock
[Yrand, in_param,out_param] = meanMCBernoulli_g_param(varargin{:});
out_param.n_hoeff = ceil(log(2/in_param.alpha)/(2*in_param.abstol^2));
out_param.n_clt = ceil(stdnorminv(1-in_param.alpha/2)/(4*in_param.abstol^2));
out_param.n = min(max(out_param.n_hoeff,out_param.n_clt),in_param.nmax);
%%  Split The Param.n into columns
nopt=min(in_param.npcmax,out_param.n);
% numbers of samples per loop step
nn=floor(out_param.n/nopt); % number of loop steps
nremain=out_param.n-nn*nopt;
% number of samples in last loop step
nloop=repmat(nopt,1,nn);
%vector of numbers of samples per loop step
if nremain>0; nloop=[nloop nremain]; nn=nn+1; end
sumY=0;
for iloop=1:nn %loops to save memory
    sumY=sumY+sum(Yrand(nloop(iloop)));
end
%%  Estimate mu
out_param.mu=sumY/out_param.n; %calculate the mean
mu=out_param.mu; %assign answer
out_param.time=toc(tstart); %elapsed time

function x = stdnorminv(p)
% this function is the inverse function of CDF of standard normal distribution
x = -sqrt(2).*erfcinv(2*p);
end

function  [Yrand, in_param,out_param] = meanMCBernoulli_g_param(varargin)
default.reltol = 1e-1;% default relative error tolerance
default.abstol  = 1e-2;% default absolute error tolerance
default.alpha = 0.01;% default uncertainty
default.npcmax = 1e6;% default n piece maximum
default.nmax = 1e8;% default n maximum
if isempty(varargin)
    help meanMCBernoulli_g
    warning('MATLAB:meanMCBernoulli_g:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand = (binornd(1,0.5,n,1)).^2.')
    p = 0.5;
    Yrand = @(n) (binornd(1,p,n,1)).^2;
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
    in_param.reltol = default.reltol;    
    in_param.alpha = default.alpha;
    in_param.npcmax = default.npcmax;
    in_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'Yrand',@isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);        
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'npcmax',default.npcmax,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end     
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'reltol',default.reltol,@isnumeric);          
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'npcmax',default.npcmax,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,Yrand,varargin{2:end})
    in_param = p.Results;
end
    if (in_param.abstol <= 0) % absolute error tolerance
        warning('MATLAB:meanMCBernoulli_g:abstolneg',...
            ['Absolute error tolerance should be greater than 0, ' ...
            'use the absolute value of the error tolerance'])
        in_param.abstol = abs(in_param.abstol);
    end
    if (in_param.reltol <= 0 || in_param.reltol >=1) % relative error tolerance
        warning('MATLAB:meanMCBernoulli_g:reltolnotin01',...
            ['Relative error tolerance should be less than 1 and bigger than 0, ' ...
            'use the default value of the error tolerance'])
        in_param.reltol = default.reltol;
    end
    if (in_param.alpha <= 0 ||in_param.alpha >= 1) % uncertainty
        warning('MATLAB:meanMCBernoulli_g:alphanotin01',...
            ['the uncertainty should be less than 1 and bigger than 0, '...
            'use the default value.'])
        in_param.alpha = default.alpha;
    end
    if (~isposint(in_param.npcmax))
        % maxinum number of scalar values of x per vector should be a positive integer
        warning('MATLAB:meanMCBernoulli_g:npcmaxnotposint',...
            ['the number of each piece of the samples should be' ...
            'a positive integer, take the absolute value and ceil.'])
        in_param.npcmax = ceil(abs(in_param.npcmax));
    end
    if (~isposint(in_param.nmax)) % sample budget should be a positive integer
        warning('MATLAB:meanMCBernoulli_g:nmaxnotposint',...
            ['the number of nmax should be a positive integer,'...
            'take the absolute value and ceil.'])
        in_param.nbudget = ceil(abs(in_param.nmax));
    end
    out_param = in_param;%let the out_param contains all the in_param
end
end