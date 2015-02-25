function [pHat,out_param]=meanMCBer_g(varargin)
%MEANMCBER_G Monte Carlo method to estimate the mean of a Bernoulli random
%variable to within a specified absolute error tolerance with guaranteed
%confidence level 1-alpha.
%
%   pHat = MEANMCBER_G(Yrand) estimates the mean of a Bernoulli random
%   variable Y to within a specified absolute error tolerance with
%   guaranteed confidence level 99%. Input Yrand is a function handle that
%   accepts a positive integer input n and returns a n x 1 vector of IID
%   instances of the Bernoulli random variable Y.
% 
%   pHat = MEANMCBER_G(Yrand,abstol,alpha,nmax) estimates the mean
%   of a Bernoulli random variable Y to within a specified absolute error
%   tolerance with guaranteed confidence level 1-alpha using all ordered
%   parsing inputs abstol, alpha and nmax.
% 
%   pHat = MEANMCBER_G(Yrand,'abstol',abstol,'alpha',alpha,'nmax',nmax)
%   estimates the mean of a Bernoulli random variable Y to within a
%   specified absolute error tolerance with guaranteed confidence level
%   1-alpha. All the field-value pairs are optional and can be supplied in
%   different order.
% 
%   [pHat, out_param] = MEANMCBER_G(Yrand,in_param) estimates the mean
%   of a Bernoulli random variable Y to within a specified absolute error
%   tolerance with the given parameters in_param and produce the estimated
%   mean pHat and output parameters out_param.
% 
%   Input Arguments
%
%     Yrand --- the function for generating IID instances of a Bernoulli
%               random variable Y whose mean we want to estimate.
%
%     pHat --- the estimated mean of Y.
%
%     in_param.abstol --- the absolute error tolerance, the default value is 1e-2.
% 
%     in_param.alpha --- the uncertainty, the default value is 1%.
% 
%     in_param.nmax --- the sample budget, the default value is 1e9.
% 
%   Output Arguments
%
%     out_param.n --- the total sample used.
%
%     out_param.time --- the time elapsed in seconds.
% 
%  Guarantee
%
% If the sample size is calculated according Hoeffding's inequality, which
% equals to ceil(log(2/out_param.alpha)/(2*out_param.abstol^2)), then the
% following inequality must be satisfied:
%
% Pr(|p-pHat| <= abstol) >= 1-alpha.
% 
% Here p is the true mean of Yrand, and pHat is the output of MEANMCBER_G
%
% Also, the cost is deterministic.
%
%   Examples
% 
%   Example 1:
%   Calculate the mean of a Bernoulli random variable with true p=1/90,
%   absolute error tolerance 1e-3 and uncertainty 0.01.
% 
%   >> in_param.abstol=1e-3; in_param.alpha = 0.01; p=1/9;Yrand=@(n) rand(n,1)<p;
%   >> pHat=meanMCBer_g(Yrand,in_param)
%   pHat = 0.1***
% 
% 
%   Example 2:
%   Using the same function as example 1, with the absolute error tolerance
%   1e-4.
% 
%   >> pHat = meanMCBer_g(Yrand,1e-3)
%   pHat = 0.11***
% 
% 
%   Example 3:
%   Using the same function as example 1, with the absolute error
%   tolerance 1e-2 and uncertainty 0.05.
% 
%   >> pHat = meanMCBer_g(Yrand,'abstol',1e-2,'alpha',0.05)
%   pHat = 0.1***
% 
% 
%   See also FUNAPPX_G, INTEGRAL_G, CUBMC_G, MEANMC_G, CUBLATTICE_G, CUBSOBOL_G
% 
%  References
%
%   [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
%   Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014.
%   arXiv:1208.4318 [math.ST]
%
%   [2] Sou-Cheng T.  Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   "GAIL: Guaranteed Automatic Integration Library (Version 2.1)"
%   [MATLAB Software], 2015. Available from
%   http://code.google.com/p/gail/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software", Journal of Open Research
%   Software, Volume 2, Number 1 (2014), e22, pp. 1-7, DOI:
%   http://dx.doi.org/10.5334/jors.bb (describes principles of Reliable
%   Reproducible Research and Supportable Scientific Software)
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. (develops practices of Reliable
%   Reproducible Research and Supportable Scientific Software)
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

tstart = tic; %start the clock

[Yrand,out_param] = meanMCBer_g_param(varargin{:});
%parameter checking and parsing
out_param.exit = 0;
out_param = nabs(out_param);%calculate the sample needed
npcmax = 1e6;
out_param.pHat = gail.evalmean(Yrand,out_param.n,npcmax);
% evaluate the mean 
pHat=out_param.pHat; % assign answer
out_param.time=toc(tstart); %elapsed time
end

function  [Yrand,out_param] = meanMCBer_g_param(varargin)
default.abstol  = 1e-2;% default absolute error tolerance
default.alpha = 1e-2;% default uncertainty
default.nmax = 1e9;% default the sample budget
if isempty(varargin) % if no value was parsed
    help meanMCBer_g % print documentation
    warning('MATLAB:meanMCBer_g:yrandnotgiven',...
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

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
  f_addParamVal= @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin
    %if there is only one input which is Yrand, use all the default parameters
    out_param.abstol = default.abstol;
    out_param.alpha = default.alpha;
    out_param.nmax = default.nmax;
else
    pHat = inputParser;
    addRequired(pHat,'Yrand',@gail.isfcn);
    if isnumeric(in2)
    %if there are multiple inputs with only numeric, they should be put in order.
        addOptional(pHat,'abstol',default.abstol,@isnumeric);
        addOptional(pHat,'alpha',default.alpha,@isnumeric);
        addOptional(pHat,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            pHat.StructExpand = true;
            pHat.KeepUnmatched = true;
        end
        f_addParamVal(pHat,'abstol',default.abstol,@isnumeric);
        f_addParamVal(pHat,'alpha',default.alpha,@isnumeric);
        f_addParamVal(pHat,'nmax',default.nmax,@isnumeric);
    end
    parse(pHat,Yrand,varargin{2:end})
    out_param = pHat.Results;
end
if (~gail.isfcn(Yrand))
    warning('MATLAB:meanMCBer_g:yrandnotfcn',...
        ['The input must be a function handle.',...
        ' Now GAIL is using default Bernoulli random variable with parameter 0.0078.'])%print warning message
    p = 2^(-7);
    Yrand = @(n) rand(n,1)<p;% use the default random variable
end
if max(size(Yrand(5)))~=5 || min(size(Yrand(5)))~=1
    % if the input is not a length n vector, print the warning message
    warning('MATLAB:meanMCBer_g:yrandnotlengthN',...
        ['Yrand should be a random variable vector of length n, '...
        'but not an integrand or a matrix. Now GAIL is using the default Yrand.'])
    p = 2^(-7);
    Yrand = @(n) rand(n,1)<p;% use the default random variable
end
if (out_param.abstol < 0) %absolute error tolerance negative
    warning('MATLAB:meanMCBer_g:abstolneg',...
        ['Absolute error tolerance should be greater than 0; ' ...
        'We will use the absolute value of the error tolerance'])
    out_param.abstol = abs(out_param.abstol);
end
if (out_param.alpha <= 0 ||out_param.alpha >= 1) 
    % uncertainty is not in (0,1)
    warning('MATLAB:meanMCBer_g:alphanotin01',...
        ['the uncertainty should be in (0,1); '...
        'We will use the default value.'])
    out_param.alpha = default.alpha;
end

if (~gail.isposge30(out_param.nmax)) 
    % sample budget should be a large positive integer
    warning('MATLAB:meanMCBer_g:nmaxnotposint',...
        ['the number of nmax should be a large positive integer; '...
        'We will take the default value 1e9.'])
    out_param.nbudget = default.nmax;
end
end

function out_param=meanMCBernoulli_g_err(out_param)
% Handles errors in meanMCBer_g to give an exit with information.
% out_param.exit = 0   success
%                  1   n exceed nmax
if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
    case 1 % nabs exceed nmax.
        warning('MATLAB:meanMCBer_g:nabsexceednmax',...
            [' To guarantee the absolute error, tried to evaluate at '...
            int2str(out_param.n) ...
            ' samples, which is more than the allowed maximum of '...
            num2str(out_param.nmax) ...
            ' samples. We will use the maximum sample budget to estimate p '...
            'without guarantee.']);
        return;
end
end

function out_param = nabs(out_param)
% this function is to calculate the sample size n needed to satisfy
% absolute error criterion
out_param.n = ceil(log(2/out_param.alpha)/(2*out_param.abstol^2));
% calculate the sample size by Hoeffding's inequality
out_param.tau = 1;
% it is one step estimation
if out_param.n > out_param.nmax % if the sample needed is bigger than nmax
    out_param.exit=1; % pass a flag
    meanMCBernoulli_g_err(out_param);% print warning message
    out_param.n = out_param.nmax;% update nabs
end
end

