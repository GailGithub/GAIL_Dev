function [tmu,out_param]=meanMC_g(varargin)
% meanMC_G Monte Carlo method to estimate the mean of a random variable.
%
%   tmu = meanMC_G(Yrand) estimates the mean, mu, of a random variable Y to
%   within a specified generalized error tolerance, tolfun :=
%   max(abstol,reltol*|mu|), i.e., |mu - tmu| <= tolfun with probability at
%   least (1 - alpha), where abstol is the absolute error tolerance, and
%   reltol is the relative error tolerance. Usually the reltol determines
%   the accuracy of the estimation, however, if |mu| is rather small, then
%   abstol determines the accuracy of the estimation. Input Yrand is a
%   function handle that accepts a positive integer input n and returns an
%   n x 1 vector of IID instances of the random variable Y.
%
%   tmu = meanMC_G(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%   estimates the mean of a random variable Y to within a specified
%   generalized error tolerance tolfun with guaranteed confidence
%   level 1-alpha using all ordered parsing inputs abstol, reltol, alpha,
%   fudge, nSig, n1, tbudget, nbudget.
%
%   tmu = meanMC_G(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha,
%   'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',tbudget,'nbudget',nbudget)
%   estimates the mean of a random variable Y to within a specified
%   generalized error tolerance tolfun with guaranteed confidence level
%   1-alpha. All the field-value pairs are optional and can be supplied in
%   different order, if a field is not supplied, the default value is used.
%
%   [tmu, out_param] = meanMC_G(Yrand,in_param) estimates the mean of a
%   random variable Y to within a specified generalized error tolerance
%   tolfun with the given parameters in_param and produce the estimated
%   mean tmu and output parameters out_param. If a field is not specified,
%   the default value is used.
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
%     in_param.abstol --- the absolute error tolerance, which should be
%     positive, default value is 1e-2.
%
%     in_param.reltol --- the relative error tolerance, which should be
%     between 0 and 1, default value is 1e-1.
%
%     in_param.alpha --- the uncertainty, which should be a small positive
%     percentage, default value is 1%.
%
%     in_param.fudge --- standard deviation inflation factor, which should
%     be larger than 1, default value is 1.2.
%
%     in_param.nSig --- initial sample size for estimating the sample
%     variance, which should be a moderately large integer bigger than or
%     equal to 30, the default value is 1e4.
%
%     in_param.n1 --- initial sample size for estimating the sample mean,
%     which should be a moderately large positive integer at least 30, the
%     default value is 1e4.
%
%     in_param.tbudget --- the time budget in seconds to do the two-stage
%     estimation, which should be positive, the default value is 100 seconds.
%
%     in_param.nbudget --- the sample budget to do the two-stage
%     estimation, which should be a large positive integer, the default
%     value is 1e9.
%
%   Output Arguments
%
%     tmu --- the estimated mean of Y.
%
%     out_param.tau --- the total number of iterations.
%
%     out_param.n --- the sample size used in each iteration.
%
%     out_param.nremain --- the remaining sample budget to estimate mu. It was
%     calculated by the sample left and time left.
%
%     out_param.ntot --- total sample used, including the sample used to
%     convert time budget to sample budget and the sample in each iteration
%     step.
%
%     out_param.hmu --- estimated mean in each iteration.
%
%     out_param.bound_err --- a reliable upper bound on error for each iteration.
%
%     out_param.var --- the sample variance.
%
%     out_param.exitflag --- the state of program when exiting.
%
%                      0   Success
%
%                      1   Not enough samples to estimate the mean
%
%     out_param.kurtmax --- the upper bound on modified kurtosis.
%
%     out_param.time --- the time elapsed in seconds.
%
%
%  Guarantee
% This algorithm attempts to calculate the mean, mu, of a random variable
% to a prescribed error tolerance, tolfun:= max(abstol,reltol*|mu|), with
% guaranteed confidence level 1-alpha. If the algorithm terminates without
% showing any warning messages and provides an answer tmu, then the follow
% inequality would be satisfied:
%
% Pr(|mu-tmu| <= tolfun) >= 1-alpha
%
% The cost of the algorithm, N_tot, is also bounded above by N_up, which is
% defined in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And
% the following inequality holds:
%
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
%
% Examples
%
% Example 1:
% If no parameters are parsed, help text will show up as follows:
% >> meanMC_g
% ***Monte Carlo method to estimate***
%
%
% Example 2:
% Calculate the mean of x^2 when x is uniformly distributed in
% [0 1], with the absolute error tolerance = 1e-3 and uncertainty 5%.
%
% >> in_param.reltol = 0; in_param.abstol = 1e-3;
% >> in_param.alpha = 0.05; Yrand=@(n) rand(n,1).^2;
% >> tmu = meanMC_g(Yrand,in_param); exactsol = 1/3;
% >> check = double(abs(exactsol-tmu) < 1e-3)
% check = 1
%
%
% Example 3:
% Calculate the mean of exp(x) when x is uniformly distributed in
% [0 1], with the absolute error tolerance 1e-3.
%
% >> tmu = meanMC_g(@(n)exp(rand(n,1)),1e-3,0); exactsol = exp(1)-1;
% >> check = double(abs(exactsol-tmu) < 1e-3)
% check = 1
%
%
% Example 4:
% Calculate the mean of cos(x) when x is uniformly distributed in
% [0 1], with the relative error tolerance 1e-2 and uncertainty 0.05.
%
% >> tmu = meanMC_g(@(n)cos(rand(n,1)),'reltol',1e-3,'abstol',1e-4,'alpha',0.01);
% >> exactsol = sin(1);
% >> check = double(abs(exactsol-tmu) < max(1e-3,1e-2*abs(exactsol)))
% check = 1
%
%
%   See also FUNAPPX_G, INTEGRAL_G, CUBMC_G, CUBSOBOL_G, CUBLATTICE_G, CUBBAYESLATTICE_G
%
%  References
%
%   [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, "Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo
%   sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
%   Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), pp. 105-128,
%   Springer-Verlag, Berlin, 2014. DOI: 10.1007/978-3-642-41095-6_5
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice and Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   [6] Fang, K.-T., & Wang, Y. (1994). Number-theoretic Methods in
%   Statistics. London, UK: CHAPMAN & HALL
%
%   [7] Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
%   Means of Random Variables, PhD Thesis, Illinois Institute of
%   Technology, 2016.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above paper and software.
%
%   Author:  Lan Jiang

tstart = tic; %start the clock
[Yrand, out_param] = meanMC_g_param(varargin{:});
out_param.Yrand = Yrand;
n1 = 2;
Yrand(n1); %let it run once to load all the data. warm up the machine.
nsofar = n1;

ntry = 10;% try several samples to get the time
tic;
Yrand(ntry);
ttry=toc;
tpern = ttry/ntry; % calculate time per sample
nsofar = nsofar+ntry; % update n so far
out_param.exitflag = 0;
if tpern<1e-7;%each sample use very very little time
    booster = 8;
    tic;Yrand(ntry*booster);ttry2 = toc;
    ntry = ntry*[1 booster];
    ttry = [ttry ttry2];% take eight times more samples to try
elseif tpern>=1e-7 && tpern<1e-5 %each sample use very little time
    booster = 6;
    tic;Yrand(ntry*booster);ttry2 = toc;
    ntry = ntry*[1 booster];
    ttry = [ttry ttry2];% take six times more samples to try
elseif tpern>=1e-5 && tpern<1e-3 %each sample use little time
    booster = 4;
    tic;Yrand(ntry*booster);ttry2 = toc;
    ntry = ntry*[1 booster];
    ttry = [ttry ttry2];% take four times more samples to try
elseif  tpern>=1e-3 && tpern<1e-1 %each sample use moderate time
    booster = 2;
    tic;Yrand(ntry*booster);ttry2 = toc;
    ntry = ntry*[1 booster];
    ttry = [ttry ttry2];% take two times more samples to try
else %each sample use lots of time, stop try
end
[tmu,out_param] =  meanmctolfun(Yrand,out_param,ntry,ttry,nsofar,tstart);
%control the order of out_param
if out_param.reltol ~= 0
    out_param = orderfields(out_param, ...
        {'Yrand','abstol','reltol','bound_err','alpha','fudge', 'tau','hmu','time',...
        'n1','nSig', 'n','nremain','nbudget','ntot','tbudget','var',...
        'kurtmax','exitflag'});
else
    out_param = orderfields(out_param, ...
        {'Yrand','abstol','reltol','bound_err','alpha','fudge', 'tau','time',...
        'n1','nSig', 'n','nremain','nbudget','ntot','tbudget','var',...
        'kurtmax','exitflag'}); % no hmu
end
end

function [tmu,out_param] =  meanmctolfun(Yrand,out_param,ntry,ttry,nsofar,tstart)
tic;
Yval = Yrand(out_param.nSig);% get samples to estimate variance
t_sig = toc;%get the time for calculating nSig function values.
nsofar = nsofar+out_param.nSig;% update the samples that have been used
out_param.nremain = gail.estsamplebudget(out_param.tbudget,...
    out_param.nbudget,[ntry out_param.nSig],nsofar,tstart,[ttry t_sig]);
%update the nremain could afford until now
out_param.var = var(Yval);% calculate the sample variance--stage 1
sig0 = sqrt(out_param.var);% standard deviation
sig0up = out_param.fudge*sig0;% upper bound on the standard deviation
alpha_sig = out_param.alpha/2;% the uncertainty for variance estimation
out_param.kurtmax = (out_param.nSig-3)/(out_param.nSig-1) ...
    + ((alpha_sig*out_param.nSig)/(1-alpha_sig))...
    *(1-1/out_param.fudge^2)^2;
%the upper bound on the modified kurtosis
npcmax = 1e6;%constant to do iteration and mean calculation4
if out_param.reltol == 0
    out_param.tau = 1;
    alphai = 1-(1-out_param.alpha)/(1-alpha_sig);
    if sig0up == 0 % if the variance is zero, just take n_sigma samples
        out_param.n = out_param.nSig;
    else
        toloversig = out_param.abstol/sig0up;
        % absolute error tolerance over sigma
        [out_param.n,out_param.bound_err] = nchebe(toloversig,alphai,out_param.kurtmax,out_param.nbudget,sig0up);
        if out_param.n > out_param.nremain;
            out_param.exitflag=1; %pass a flag
            meanMC_g_err(out_param); % print warning message
            out_param.n = out_param.nremain;% update n
        end
    end
    tmu = gail.evalmean(Yrand,out_param.n,npcmax);%evaluate the mean
    nsofar = nsofar+out_param.n;%total sample used
else
    alphai = (out_param.alpha-alpha_sig)/(2*(1-alpha_sig));
    %uncertainty to do iteration
    eps1 = ncbinv(out_param.n1,alphai,out_param.kurtmax);
    %tolerance for initial estimation
    out_param.bound_err(1) = sig0up*eps1;
    %the width of initial confidence interval for the mean
    i=1;
    out_param.n(i) = out_param.n1;% initial sample size to do iteration
    while true
        out_param.tau = i;%step of the iteration
        if out_param.n(i) > out_param.nremain;
            % if the sample size used for initial estimation is
            % larger than nremain, print warning message and use nremain
            out_param.exitflag=1; %pass a flag
            meanMC_g_err(out_param); % print warning message
            out_param.n(i) = out_param.nremain;% update n
            tmu = gail.evalmean(Yrand,out_param.n(i),npcmax);%evaluate the mean
            nsofar = nsofar+out_param.n(i);%total sample used
            break;
        end
        out_param.hmu(i) = gail.evalmean(Yrand,out_param.n(i),npcmax);%evaluate mean
        nsofar = nsofar+out_param.n(i);
        out_param.nremain = out_param.nremain-out_param.n(i);%update n so far and nremain
        errtype = 'max';
        % error type, see the function 'tolfun' at Algoithms/+gail/ directory
        % for more info
        theta  = 0;% relative error case
        deltaplus = (gail.tolfun(out_param.abstol,out_param.reltol,...
            theta,out_param.hmu(i) - out_param.bound_err(i),errtype)...
            +gail.tolfun(out_param.abstol,out_param.reltol,...
            theta,out_param.hmu(i) + out_param.bound_err(i),errtype))/2;
        % a combination of tolfun, which used to decide stopping time
        if deltaplus >= out_param.bound_err(i) % stopping criterion
            deltaminus= (gail.tolfun(out_param.abstol,out_param.reltol,...
                theta,out_param.hmu(i) - out_param.bound_err(i),errtype)...
                -gail.tolfun(out_param.abstol,out_param.reltol,...
                theta,out_param.hmu(i) + out_param.bound_err(i),errtype))/2;
            % the other combination of tolfun, which adjust the hmu a bit
            tmu = out_param.hmu(i)+deltaminus;
            break;
        else
            out_param.bound_err(i+1) = min(out_param.bound_err(i)/2,max(out_param.abstol,...
                0.95*out_param.reltol*abs(out_param.hmu(i))));
            i=i+1;
        end
        toloversig = out_param.bound_err(i)/sig0up;%next tolerance over sigma
        alphai = (out_param.alpha-alpha_sig)/(1-alpha_sig)*2.^(-i);
        %update the next uncertainty
        out_param.n(i) = nchebe(toloversig,alphai,out_param.kurtmax,out_param.nbudget,sig0up);
    end
end
%get the next sample size needed
out_param.ntot = nsofar;%total sample size used
out_param.time=toc(tstart); %elapsed time
end

function [ncb, err] = nchebe(toloversig,alpha,kurtmax,nbudget,sig0up)
%this function uses Chebyshev and Berry-Esseen Inequality to calculate the
%sample size needed
ncheb = ceil(1/(toloversig^2*alpha));%sample size by Chebyshev's Inequality
A=18.1139;
A1=0.3328;
A2=0.429; % three constants in Berry-Esseen inequality
M3upper = kurtmax^(3/4);
%the upper bound on the third moment by Jensen's inequality
BEfun2=@(logsqrtn)gail.stdnormcdf(-exp(logsqrtn).*toloversig)...
    +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
    A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))-alpha/2;
% Berry-Esseen function, whose solution is the sample size needed
logsqrtnCLT=log(gail.stdnorminv(1-alpha/2)/toloversig);%sample size by CLT
nbe=ceil(exp(2*fzero(BEfun2,logsqrtnCLT)));
%calculate Berry-Esseen n by fzero function
ncb = min(min(ncheb,nbe),nbudget);%take the min of two sample sizes.
if(nargout>1)
    logsqrtn = log(sqrt(ncb));
    BEfun3=@(toloversig)gail.stdnormcdf(-exp(logsqrtn).*toloversig)...
        +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
        A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))-alpha/2;
    err = fzero(BEfun3,toloversig) * sig0up;
end
end

function eps = ncbinv(n1,alpha1,kurtmax)
%This function calculate the reliable upper bound on error when given
%Chebyshev and Berry-Esseen inequality and sample size n.
NCheb_inv = 1/sqrt(n1*alpha1);
% use Chebyshev inequality
A=18.1139;
A1=0.3328;
A2=0.429; % three constants in Berry-Esseen inequality
M3upper=kurtmax^(3/4);
%using Jensen's inequality to bound the third moment
BEfun=@(logsqrtb)gail.stdnormcdf(n1.*logsqrtb)...
    +min(A1*(M3upper+A2), ...
    A*M3upper./(1+(sqrt(n1).*logsqrtb).^3))/sqrt(n1)...
    - alpha1/2;
% Berry-Esseen inequality
logsqrtb_clt=log(sqrt(gail.stdnorminv(1-alpha1/2)/sqrt(n1)));
%use CLT to get tolerance
NBE_inv = exp(2*fzero(BEfun,logsqrtb_clt));
%use fzero to get Berry-Esseen tolerance
eps = min(NCheb_inv,NBE_inv);
%take the min of Chebyshev and Berry Esseen tolerance
end

function  [Yrand,out_param] = meanMC_g_param(varargin)

default.abstol  = 1e-2;% default absolute error tolerance
default.reltol = 1e-1;% default relative error tolerance
default.nSig = 1e4;% default initial sample size nSig for variance estimation
default.n1 = 1e4; % default initial sample size n1 for mean estimation
default.alpha = 0.01;% default uncertainty
default.fudge = 1.2;% default fudge factor
default.tbudget = 100;% default time budget
default.nbudget = 1e9; % default sample budget
if isempty(varargin)
    help meanMC_g
    warning('GAIL:meanMC_g:yrandnotgiven',...
        'Yrand must be specified. Now GAIL is using Yrand =@(n) rand(n,1).^2.')
    Yrand = @(n) rand(n,1).^2;
    out_param.Yrand = Yrand;
    %if no values are parsed, print warning message and use the default
    %random variable
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
    %if there is only input which is Yrand, use all the default parameters
    out_param.abstol = default.abstol;% default absolute error tolerance
    out_param.reltol = default.reltol; % default relative error tolerance
    out_param.alpha = default.alpha;% default uncertainty
    out_param.fudge = default.fudge;% default standard deviation inflation factor
    out_param.nSig = default.nSig;% default the sample size to estimate the variance
    out_param.n1 = default.n1;% default the initial sample size to estimate the mean
    out_param.tbudget = default.tbudget;% default time budget
    out_param.nbudget = default.nbudget;% default sample budget
else
    p = inputParser;
    addRequired(p,'Yrand',@gail.isfcn);
    if isnumeric(in2)
        %if there are multiple inputs with only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'nSig',default.nSig,@isnumeric);
        addOptional(p,'n1',default.n1,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'nSig',default.nSig,@isnumeric);
        addParamValue(p,'n1',default.n1,@isnumeric);
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
    end
    parse(p,Yrand,varargin{2:end})
    out_param = p.Results;
end
if (~gail.isfcn(Yrand))
    warning('GAIL:meanMC_g:yrandnotfcn',...
        ['Yrand must be a function handle.'...
        ' Now GAIL is using default Yrand =@(n) rand(n,1).^2 .'])
    %print warning message
    Yrand = @(n) rand(n,1).^2;
end
if max(size(Yrand(5)))~=5 || min(size(Yrand(5)))~=1
    % if the input is not a length n vector, print the warning message
    warning('GAIL:meanMC_g:yrandnotlengthN',...
        ['Yrand should be a random variable vector of length n, '...
        'but not an integrand or a matrix.'...
        ' Now GAIL is using the default Yrand =@(n) rand(n,1).^2.'])
    Yrand = @(n) rand(n,1).^2;
end

if (out_param.abstol < 0)
    %absolute error tolerance should be positive
    warning('GAIL:meanMC_g:abstolneg',...
        ['Absolute error tolerance should be greater than 0; ' ...
        'We will take the aosolute value of the absolute error tolerance provided.'])
    out_param.abstol = abs(out_param.abstol);
end
if (out_param.reltol < 0 || out_param.reltol > 1)
    % relative error tolerance should be in [0,1]
    warning('GAIL:meanMC_g:reltolneg',...
        ['Relative error tolerance should be between 0 and 1; ' ...
        'We will use the default value of the relative error tolerance 1e-1.'])
    out_param.reltol = default.reltol;
end
if (out_param.alpha <= 0 ||out_param.alpha >= 1)
    %uncertainty should be in (0,1)
    warning('GAIL:meanMC_g:alphanotin01',...
        ['The uncertainty should be between 0 and 1; '...
        'We will use the default value 0.01.'])
    out_param.alpha = default.alpha;
end
if (out_param.fudge <= 1)
    %standard deviation inflation factor should be bigger than 1
    warning('GAIL:meanMC_g:fudgelessthan1',...
        ['The fudge factor should be larger than 1; '...
        'We will use the default value 1.2.'])
    out_param.fudge = default.fudge;
end
if (~gail.isposge30(out_param.nSig))
    %initial sample size should be a positive integer at least 30
    warning('GAIL:meanMC_g:nsignotposint',...
        ['The number nSig should a positive integer at least 30; '...
        'We will use the default value 1e4.'])
    out_param.nSig = default.nSig;
end
if (~gail.isposge30(out_param.n1))
    %initial sample size should be a posotive integer at least 30
    warning('GAIL:meanMC_g:n1notposint',...
        ['The number n1 should a positive integer at least 30; '...
        'We will use the default value 1e4.'])
    out_param.n1 = default.n1;
end
if (out_param.tbudget < 0)
    %time budget in seconds should be positive
    warning('GAIL:meanMC_g:timebudgetneg',...
        ['Time budget in seconds should be positive; '...
        'We will take the absolute value of time budget provided'])
    out_param.tbudget = abs(out_param.tbudget);
end

if (out_param.tbudget == 0)
    %time budget in seconds should be positive
    warning('GAIL:meanMC_g:timebudget0',...
        ['Time budget in seconds should be positive rather than zero; '...
        'We will take the default value of time budget 100 seconds'])
    out_param.tbudget = abs(out_param.tbudget);
end
if (~gail.isposge30(out_param.nbudget))
    %sample budget should be a large positive integer
    warning('GAIL:meanMC_g:nbudgetnotposint',...
        ['The number of sample budget should be a large positive integer; '...
        'We will use the default value 1e9.'])
    out_param.nbudget =default.nbudget;
end
out_param.bound_err = inf; %create this field
out_param.exitflag = 1;
%pass the signal indicating the parameters have been checked
end

function out_param = meanMC_g_err(out_param)
% Handles errors in meanMC_g and meanMC_g_param to give an exit with
%  information.
%            out_param.exitflag = 0   success
%                             1   too many samples required

if ~isfield(out_param,'exitflag'); return; end
if out_param.exitflag==0; return; end
switch out_param.exitflag
    case 1 % not enough samples to estimate the mean.
        nexceed = out_param.n(out_param.tau);
        warning('GAIL:meanMC_g:maxreached',...
            [' In order to achieve the guaranteed accuracy, at step '...
            int2str(out_param.tau) ' , tried to evaluate at ' int2str(nexceed) ...
            ' samples, which is more than the remaining '...
            int2str(out_param.nremain) ...
            ' samples. We will use all the samples left to estimate the mean without guarantee.']);
        return
end
end
