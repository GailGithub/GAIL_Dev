function [pHat,out_param]=meanMCBernoulli_g(varargin)
%MEANMCBERNOULLI_G Monte Carlo method to estimate the mean of a Bernoulli
%random variable to within a specified error tolerance with guaranteed
%confidence level 1-alpha.
%
%   pHat = MEANMCBERNOULLI_G(Yrand) estimates the mean of a Bernoulli random
%   variable Y to within a specified error tolerance with guaranteed
%   confidence level 99%. Input Yrand is a function handle that accepts a
%   positive integer input n and returns a n x 1 vector of IID instances
%   of the Bernoulli random variable Y.
%
%   pHat = MEANMCBERNOULLI_G(Yrand,abstol,reltol,errtype,alpha,nmax) estimates
%   the mean of a Bernoulli random variable Y to within a specified error
%   tolerance with guaranteed confidence level 1-alpha using all ordered
%   parsing inputs abstol, reltol, errtype, alpha and nmax.
%
%   pHat = MEANMCBERNOULLI_G(Yrand,'abstol',abstol,'reltol',reltol,'errtype',errtype,'alpha',alpha,'nmax',nmax)
%   estimates the mean of a Bernoulli random variable Y to within a
%   specified error tolerance with guaranteed confidence level 1-alpha. All
%   the field-value pairs are optional and can be supplied in different order.
%
%   [pHat, out_param] = MEANMCBERNOULLI_G(Yrand,in_param) estimates the mean
%   of a Bernoulli random variable Y to within a specified error tolerance
%   with the given parameters in_param and produce the estimated mean pHat
%   and output parameters out_param.
%
%   Input Arguments
%
%     Yrand --- the function for generating IID instances of a Bernoulli
%               random variable Y whose mean we want to estimate.
%
%     pHat --- the estimated mean of Y.
%
%     in_param.abstol --- the absolute error tolerance, default value is 1e-2.
%
%     in_param.reltol --- the relative error tolerance, default value is 1e-1.
%
%     in_param.errtype --- the error type, default value is 'abs'.
%
%                          'abs'--- absolute error criterion
%
%                          'rel'--- relative error criterion
%
%                          'either'---absolute OR relative criterion
%
%     in_param.alpha --- the uncertainty, default value is 1%.
%
%     in_param.nmax --- the sample budget, default value is 1e9.
%
%   Output Arguments
%
%     out_param.nabs --- sample size needed to satisfy absolute error
%                        tolerance.
%
%     out_param.nrel --- sample size needed to satisfy relative error
%                        tolerance.
%
%     out_param.n --- the total sample used.
%
%     out_param.tau --- the iteration step.
%
%     out_param.time --- the time elapsed in seconds.
%
%  Guarantee
%
% Case 1: errtype = 'abs'
%
% If the sample size is calculated according Hoeffding's inequality, which
% equals to ceil(log(2/out_param.alpha)/(2*out_param.abstol^2)), then the
% following inequality must be satisfied:
%
% Pr(|p-pHat| <= abstol) >= 1-alpha.
%
% Here p is the true mean of Yrand, and pHat is the output of
% MEANMCBERNOULLI_G with errtype = 'abs'
%
% Also, the cost is deterministic and bounded.
%
% Case 2: errtype = 'rel'
%
% If the algorithm terminated without any warning messages, the estimated
% mean pHat would satisfy the following inequality:
%
% Pr(|p-pHat| <= abstol*p) >= 1-alpha.
%
% Here p is the true mean of Y, and pHat is the output of MEANMCBERNOULLI_G
% with errtype = 'rel'.
%
% Additionally, the cost of the algorithm would be bounded by N_up, which is
% defined in terms of the true mean p, uncertainty alpha and relative
% tolerance reltol. For details, please refer to the paper.
%
%   Examples
%
%   Example 1:
%   Calculate the mean of a Bernoulli random variable with true p=1/90,
%   absolute error tolerance 1e-3 and uncertainty 0.01.
%
%   >> in_param.abstol=1e-3; in_param.alpha = 0.01; p=1/9;Yrand=@(n) rand(n,1)<p;
%   >> pHat=meanMCBernoulli_g(Yrand,in_param)
%   pHat = 0.1***
%
%
%   Example 2:
%   Using the same function as example 1, with the relative error tolerance
%   1e-2.
%
%   >> pHat = meanMCBernoulli_g(Yrand,0,1e-2,'rel')
%   pHat = 0.1***
%
%
%   Example 3:
%   Using the same function as example 1, with the relative error
%   tolerance 1e-2 and uncertainty 0.05.
%
%   >> pHat = meanMCBernoulli_g(Yrand,'errtype','rel','reltol',1e-2,'alpha',0.05)
%   pHat = 0.11***
%
%
%   See also FUNAPPX_G, INTEGRAL_G, CUBMC_G, MEANMC_G, CUBLATTICE_G,
%   CUBSOBOL_G
%
%  References
%
%   [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
%   Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014.
%   arXiv:1208.4318 [math.ST]
%
%   [2] Sou-Cheng T.  Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   "GAIL: Guaranteed Automatic Integration Library (Version 2.0)"
%   [MATLAB Software], 2014. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above paper and software.
%

tstart = tic; %start the clock

[Yrand,out_param] = meanMCBernoulli_g_param(varargin{:});
%parameter checking and parsing
nsofar = 0;%sample used so far
out_param.exit = 0;
if strcmpi(out_param.errtype,'abs') % absolute error type
    out_param = nabs(out_param,tstart);%calculate the sample needed
    out_param.n = out_param.nabs;%update n
end
if strcmpi(out_param.errtype,'rel')% relative error type
    [nsofar,out_param] = nrel(out_param,Yrand,tstart);
    %calculate the sample needed
    out_param.n = out_param.nrel;%update n
end

if strcmpi(out_param.errtype,'either')
    %either absolute or relative error tolerance, which satisfies first
    out_param = nabs(out_param,tstart);
    [nsofar,out_param] = nrel(out_param,Yrand,tstart);
    out_param.n = min(out_param.nabs,out_param.nrel);
end
npcmax = 1e6;
out_param.pHat = gail.evalmean(Yrand,out_param.n,npcmax);
% evaluate the mean
pHat=out_param.pHat; % assign answer
out_param.n = out_param.n+nsofar; % update total sample used
out_param.time=toc(tstart); %elapsed time
end

function  [Yrand,out_param] = meanMCBernoulli_g_param(varargin)
default.reltol = 1e-1;% default relative error tolerance
default.abstol  = 1e-2;% default absolute error tolerance
default.errtype = 'abs';% default error type
default.alpha = 1e-2;% default uncertainty
default.nmax = 1e9;% default the sample budget
if isempty(varargin) % if no parsing value
    help meanMCBernoulli_g % print documentation
    warning('MATLAB:meanMCBernoulli_g:yrandnotgiven',...
        ['Yrand must be specified. Now GAIL is using Bernoulli random ', ...
    'variable with parameter 0.0078.'])%print warning message
    p = 2^(-7);
    Yrand = @(n) rand(n,1)<p;% use the default random variable
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
    %if there is only one input which is Yrand, use all the default parameters
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.errtype = default.errtype;
    out_param.alpha = default.alpha;
    out_param.nmax = default.nmax;
else
    pHat = inputParser;
    addRequired(pHat,'Yrand',@gail.isfcn);
    if isnumeric(in2)
    %if there are multiple inputs with only numeric, they should be put in order.
        addOptional(pHat,'abstol',default.abstol,@isnumeric);
        addOptional(pHat,'reltol',default.reltol,@isnumeric);
        addOptional(pHat,'errtype',default.errtype,@ischar);
        addOptional(pHat,'alpha',default.alpha,@isnumeric);
        addOptional(pHat,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            pHat.StructExpand = true;
            pHat.KeepUnmatched = true;
        end
        addParamValue(pHat,'abstol',default.abstol,@isnumeric);
        addParamValue(pHat,'reltol',default.reltol,@isnumeric);
        addParamValue(pHat,'errtype',default.errtype,@ischar);
        addParamValue(pHat,'alpha',default.alpha,@isnumeric);
        addParamValue(pHat,'nmax',default.nmax,@isnumeric);
    end
    parse(pHat,Yrand,varargin{2:end})
    out_param = pHat.Results;
end
if (out_param.abstol < 0) %absolute error tolerance negative
    warning('MATLAB:meanMCBernoulli_g:abstolneg',...
        ['Absolute error tolerance should be greater than 0; ' ...
        'We will use the absolute value of the error tolerance'])
    out_param.abstol = abs(out_param.abstol);
end
if (out_param.reltol < 0 || out_param.reltol >1)
    % relative error tolerance is not in [0,1]
    warning('MATLAB:meanMCBernoulli_g:reltolnotin01',...
        ['Relative error tolerance should be between 0 and 1; ' ...
        'We will use the default value of the error tolerance'])
    out_param.reltol = default.reltol;
end
if (out_param.alpha <= 0 ||out_param.alpha >= 1)
    % uncertainty is not in (0,1)
    warning('MATLAB:meanMCBernoulli_g:alphanotin01',...
        ['the uncertainty should be in (0,1); '...
        'We will use the default value.'])
    out_param.alpha = default.alpha;
end

if (~gail.isposint(out_param.nmax))
    % sample budget should be a positive integer
    warning('MATLAB:meanMCBernoulli_g:nmaxnotposint',...
        ['the number of nmax should be a positive integer; '...
        'We will take the absolute value and ceil.'])
    out_param.nbudget = ceil(abs(out_param.nmax));
end
end

function [out_param,pHat]=meanMCBernoulli_g_err(out_param,tstart)
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
            [' To guarantee the absolute error, tried to evaluate at '...
            int2str(out_param.nabs) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. We will use the maximum sample budget to estimate p '...
            'without guarantee.']);
        return;
    case 2 % ni exceed nmax.
        warning('MATLAB:meanMCBernoulli_g:niexceednmax',...
            [' To guarantee the lower bound of p, tried to evaluate at '...
            int2str(out_param.ni) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. We will use the maximum sample budget to estimate p '...
            'without guarantee.']);
        return;
    case 3 % nrel exceed nmax.
        warning('MATLAB:meanMCBernoulli_g:nrelexceednmax',...
            [' To guarantee the relative error, tried to evaluate at '...
            int2str(out_param.nrel) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. We will use the maximum sample budget to estimate p ',...
            'without guarantee.']);
        return;
end
out_param.pHat=NaN;
pHat=out_param.pHat;
if nargin>1; out_param.time=toc(tstart); end
end

function out_param = nabs(out_param,tstart)
% this function is to calculate the sample size n needed to satisfy
% absolute error criterion
out_param.nabs = ceil(log(2/out_param.alpha)/(2*out_param.abstol^2));
% calculate the sample size by Hoeffding's inequality
out_param.tau = 1;
% it is one step estimation
if out_param.nabs > out_param.nmax % if the sample needed is bigger than nmax
    out_param.exit=1; % pass a flag
    meanMCBernoulli_g_err(out_param,tstart);% print warning message
    out_param.nabs = out_param.nmax;% update nabs
end
end

function [nsofar,out_param] = nrel(out_param,Yrand,tstart)
i = 1;% initial iteration step
nsofar = 0;% sample used so far
while 1
    alphap = 1e-3;% the uncertainty for the last step to evaluate the mean
    a=2;% parameter to get uncertainty in each step
    npcmax = 1e6;
    alphai = 1-(1-out_param.alpha+alphap)^((a-1)*a^-i);
    %uncertainty in each step
    toli = out_param.reltol*2^-i;%tolerance in each step
    ni = ceil(-log(alphai)/(2*toli^2));
    %sample size obtained by one side Hoeffding's inequality
    if ni > out_param.nmax-nsofar;% if ni is bigger than sample left
        out_param.exit=2; % pass a flag
        meanMCBernoulli_g_err(out_param,tstart); % print warning message
        out_param.nrel = out_param.nmax- nsofar;%update nrel using all sample left
        break
    end
    meanP = gail.evalmean(Yrand,ni,npcmax);%evaluate mean
    nsofar = nsofar+ni;%update n used so far
    c = max(meanP-toli,0); % parameter to determine the stopping time
    delta = 1/2; % constant to determine stopping time
    if c > (meanP+toli)*delta % stopping criterion
        out_param.tau = i;%stopping time tau
        out_param.nrel = ceil(-log(alphap/2)/(2*(c*out_param.reltol)^2));
        % calculate nrel by Hoeffding's inequality
        if  out_param.nrel > out_param.nmax-nsofar
            %if nrel is bigger than sample left
            out_param.exit=3; % pass a flag
            meanMCBernoulli_g_err(out_param,tstart);% print warning message
            out_param.nrel = out_param.nmax- nsofar;%update nrel
            break
        end
        break
    else
        i=i+1;%go to the next step
    end
end
end
