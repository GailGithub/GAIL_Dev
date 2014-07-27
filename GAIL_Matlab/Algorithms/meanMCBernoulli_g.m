function [p,out_param]=meanMCBernoulli_g(varargin)
% MEANMCBERNOULLI_G Monte Carlo method to estimate the mean of a Bernoullii
% random variable to within a specific error tolerance with guaranteed
% confidence level 1-alpha.
%
% p = MEANMCBERNOULLI_G(Yrand) estimates the mean of a Bernoulli random
% variable Y to within a specified error tolerance with guaranteed
% confidence level 99%. Input Yrand is a function handle that accepts a
% positive integer input n and returns an n x 1 vector of IID instances
% of the Bernoulli random variable Y.
% 
% p = MEANMCBERNOULLI_G(Yrand,abstol,reltol,index,alpha,nmax)
% estimates the mean of a Bernoulli random variable Y to within an error
% tolerance with guaranteed confidence level 1-alpha using all ordered
% parsing inputs abstol, reltol, index, alpha and nmax.
% 
% p =
% MEANMCBERNOULLI_G(Yrand,'abstol',abstol,'reltol',reltol,'index',index,
% 'alpha',alpha,'nmax',nmax) estimates the mean of a Bernoulli random
% variable Y to within a specified error tolerance with guaranteed
% confidence level 1-alpha. All the field-value pairs are optional and
% can be supplied in different order.
% 
% p = MEANMCBERNOULLI_G(Yrand,in_param) estimates the mean of a
% Bernoulli random variable Y to within a specified error
% tolerance in_param.abstol with guaranteed uncertainty within
% in_param.alpha. If a field is not specified, the default value is used.
% 
% [p, out_param] = MEANMCBERNOULLI_G(Yrand,in_param) estimates the mean
% of a Bernoulli random variable Y to within a specified absolute error
% tolerance with the given parameters in_param and output parameters
% out_param.
% 
% 
% Yrand --- the function for generating IID instances of a Bernoulli random
% variable Y whose mean we want to estimate.
% 
% p --- the estimated mean of Y.
% 
% in_param.abstol --- the absolute error tolerance, default value is 1e-2.
% 
% in_param.reltol --- the relative error tolerance, default value is 1e-1.
% 
% in_param.index --- the error tolerance criterion, default value is
% 'abs'.
% 
% in_param.alpha --- the uncertainty, default value is 1%.
% 
% in_param.nmax --- the sample budget, default value is 1e8.
% 
% out_param.n --- the total sample size used.
% 
% out_param.time --- the time elapsed.
% 
% Examples
% 
% Example 1:
% Calculate the mean of a bernoulli random variable with true p=0.55,with
% error tolerance 1e-3 and uncertainty 0.01.
% 
% >> in_param.abstol=1e-3; in_param.alpha = 0.01; p=1/90;Yrand=@(n) rand(n,1)<p;
% >> p=meanMCBernoulli_g(Yrand,in_param)
% p = 0.01***
% 
% 
% Example 2:
% Using the same function as example 1, with the absolute error tolerance
% 1e-4.
% 
% >> p=meanMCBernoulli_g(Yrand,1e-3,1e-2,'abs')
% p = 0.01***
% 
% 
% Example 3:
% Using the sample function as example 1, with the relative error
% tolerance 1e-2 and uncertainty 0.005.
% 
% >> p=meanMCBernoulli_g(Yrand,'index','rel','reltol',1e-1,'alpha',0.05)
% p = 0.011***
% 
% 
% See also FUNAPPX_G, INTEGRAL_G, CUBMC_G, MEANMC_G
% 
% Reference
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
% conservative fixed width confidence intervals via Monte Carlo sampling,
% Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G.
% W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to
% appear, arXiv:1208.4318 [math.ST]
% 
tstart = tic; %start the clock
[Yrand,out_param] = meanMCBernoulli_g_param(varargin{:});
out_param.npcmax = 1e6;
nsofar = 0;
if strcmpi(out_param.index,'abs')
    out_param = NbyAbs(out_param,tstart);
    out_param.n = out_param.nabs;
end
if strcmpi(out_param.index,'rel')
    [nsofar,out_param] = NbyRel(out_param,Yrand,tstart);
    out_param.n = out_param.nrel;
end

if strcmpi(out_param.index,'both')
    out_param = NbyAbs(out_param,tstart);
    [nsofar,out_param] = NbyRel(out_param,Yrand,tstart);
    out_param.n = max(out_param.nabs,out_param.nrel);
end

out_param.p = SplitColumnMean(Yrand,out_param.n,out_param.npcmax);
p=out_param.p; %assign answer
out_param.n = out_param.n+nsofar;
out_param.time=toc(tstart); %elapsed time
end

function x = stdnorminv(p)
% this function is the inverse function of CDF of standard normal distribution
x = -sqrt(2).*erfcinv(2*p);
end

function p =SplitColumnMean(RV,n,npcmax)
%%  Split The Param.n into columns
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
%%  Estimate p
p=sumY/n; %calculate the mean
end

function  [Yrand,out_param] = meanMCBernoulli_g_param(varargin)

default.reltol = 1e-1;% default relative error tolerance
default.abstol  = 1e-2;% default absolute error tolerance
default.index = 'abs';
default.alpha = 1e-2;% default uncertainty
default.nmax = 1e8;% default n maximum
if isempty(varargin)
    help meanMCBernoulli_g
    warning('MATLAB:meanMCBernoulli_g:yrandnotgiven',...
        ['Yrand must be specified. Now GAIL is using Bernoulli random ', ...
    'variable with parameter 0.0078.'])
    p = 2^(-7);
    Yrand = @(n) rand(n,1)<p;
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
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.index = default.index;
    out_param.alpha = default.alpha;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'Yrand',@gail.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'index',default.index,@ischar);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'index',default.index,@ischar);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,Yrand,varargin{2:end})
    out_param = p.Results;
end
if (out_param.abstol <= 0) % absolute error tolerance
    warning('MATLAB:meanMCBernoulli_g:abstolneg',...
        ['Absolute error tolerance should be greater than 0, ' ...
        'use the absolute value of the error tolerance'])
    out_param.abstol = abs(out_param.abstol);
end
if (out_param.reltol <= 0 || out_param.reltol >=1) % relative error tolerance
    warning('MATLAB:meanMCBernoulli_g:reltolnotin01',...
        ['Relative error tolerance should be less than 1 and bigger than 0, ' ...
        'use the default value of the error tolerance'])
    out_param.reltol = default.reltol;
end
if (out_param.alpha <= 0 ||out_param.alpha >= 1) % uncertainty
    warning('MATLAB:meanMCBernoulli_g:alphanotin01',...
        ['the uncertainty should be less than 1 and bigger than 0, '...
        'use the default value.'])
    out_param.alpha = default.alpha;
end

if (~gail.isposint(out_param.nmax)) % sample budget should be a positive integer
    warning('MATLAB:meanMCBernoulli_g:nmaxnotposint',...
        ['the number of nmax should be a positive integer,'...
        'take the absolute value and ceil.'])
    out_param.nbudget = ceil(abs(out_param.nmax));
end
end

function [out_param,p]=meanMCBernoulli_g_err(out_param,tstart)
% Handles errors in meanMCBernoulli_g to give an exit with information.
% out_param.exit = 0   success
%                  1   nabs exceed nmax
%                  2   ni exceed nmax
%                  3   nrel exceed nmax
if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
    case 1 % nabs exceed nmax.
        warning('MATLAB:meanMCBernoulli_g:nabsexceednmax',...
            [' To guarantee the absolute error, tried to evalute at '...
            int2str(out_param.nabs) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. Just use the maximum sample budget to estimate p,'...
            'without guarantee.']);
        return;
    case 2 % ni exceed nmax.
        warning('MATLAB:meanMCBernoulli_g:niexceednmax',...
            [' To guarantee the lower bound of p, tried to evalute at '...
            int2str(out_param.ni) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. Just use the maximum sample budget to estimate p,'...
            'without guarantee.']);
        return;
    case 3 % nrel exceed nmax.
        warning('MATLAB:meanMCBernoulli_g:nrelexceednmax',...
            [' To guarantee the relative error, tried to evalute at '...
            int2str(out_param.nrel) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. Just use the maximum sample budget to estimate p,',...
            'without guarantee.']);
        return;        
end
out_param.p=NaN;
p=out_param.p;
if nargin>1; out_param.time=toc(tstart); end
end

function out_param = NbyAbs(out_param,tstart)
out_param.n_clt = ceil(stdnorminv(1-out_param.alpha/2)/(4*out_param.abstol^2));
out_param.n_hoeff = ceil(log(2/out_param.alpha)/(2*out_param.abstol^2));
out_param.nabs = max(out_param.n_hoeff,out_param.n_clt);
if out_param.nabs > out_param.nmax  
    out_param.exit=1; % exit the loop
    meanMCBernoulli_g_err(out_param,tstart);
    out_param.nabs = out_param.nmax;
end
end

function [nsofar,out_param] = NbyRel(out_param,Yrand,tstart)
i = 1;
nsofar = 0;
while 1
    out_param.alphap = 1e-3;
    a=2;
    out_param.alphai = 1-(1-out_param.alpha+out_param.alphap)^((a-1)*a^-i);
    out_param.toli = out_param.reltol*2^-i;
    out_param.ni = ceil(-log(out_param.alphai)/(2*out_param.toli^2));
    if out_param.ni > out_param.nmax- nsofar;
        out_param.exit=2; % exit the loop
        meanMCBernoulli_g_err(out_param,tstart); % print warning message
        %out_param.ni = out_param.nmax - nsofar; % update n_p

        out_param.nrel = out_param.nmax- nsofar;
        break;
    end
    meanY = SplitColumnMean(Yrand,out_param.ni,out_param.npcmax);
    nsofar = nsofar+out_param.ni;
    c = max(meanY-out_param.toli,0);
    delta = 1/2;
    if c > (meanY+out_param.toli)*delta
        out_param.tau = i;
        out_param.nrel = ceil(-log(out_param.alphap/2)/(2*(c*out_param.reltol)^2));
        if  out_param.nrel > out_param.nmax-nsofar
            out_param.exit=3; % exit the loop
            meanMCBernoulli_g_err(out_param,tstart);
            out_param.nrel = out_param.nmax- nsofar;
            break;
        else
            break;
        end
    else
        i=i+1;
    end
end
end