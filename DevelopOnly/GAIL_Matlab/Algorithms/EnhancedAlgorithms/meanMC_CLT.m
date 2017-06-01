function [hmu,out_param]=meanMC5_CLT(YYrand,absTol,relTol,alpha,nSig,inflate)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   tmu = MEANMC_CLT(Yrand,absTol,relTol,alpha,nSig,inflate) estimates the
%   mean, mu, of a random variable Y to within a specified error tolerance,
%   i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
%   1-alpha, where abstol is the absolute error tolerance.  The default
%   values are abstol=1e-2 and alpha=1%. Input Yrand is a function handle
%   that accepts a positive integer input n and returns an n x 1 vector of
%   IID instances of the random variable Y.
%
%   Input Arguments
%
%     Yrand --- the function or structure for generating n IID instances of a random
%     variable Y whose mean we want to estimate. Y is often defined as a
%     function of some random variable X with a simple distribution. The
%     input of Yrand should be the number of random variables n, the output
%     of Yrand should be n function values. For example, if Y = X.^2 where X
%     is a standard uniform random variable, then one may define Yrand =
%     @(n) rand(n,1).^2.
%     
%
%     absTol --- the absolute error tolerance, which should be
%     non-negative --- default = 1e-2
%
%     relTol --- the relative error tolerance, which should be
%     non-negative and no greater than 1 --- default = 0
%
%     alpha --- the uncertainty, which should be a small positive
%     percentage --- default = 1%
%
%     nSig --- the number of samples used to compute the sample variance
%     --- default = 1000
%
%     inflate --- the standard deviation inflation factor --- default = 1.2
%
%   Output Arguments
%
%     hmu --- the estimated mean of Y.
%
%     out_param.ntot --- total sample used.
%
%     out_param.var --- the sample variance.
%
%     out_param.time --- the time elapsed in seconds.
%

% This is a heuristic algorithm based on a Central Limit Theorem
% approximation
if nargin < 6
    inflate = 1.2; %standard deviation inflation factor
    if nargin < 5
        nSig = 1e3; %number of samples to estimate variance
        if nargin < 4
            alpha = 0.01; %uncertainty
            if nargin < 3
                relTol = 0.01; %relative error tolerance
                if nargin < 2
                    absTol = 1e-2; %absolute error tolerance                   
                    if nargin < 1
                        YYrand = @(n) rand(n,1); %random number generator
                       
                    end
                end
            end
        end
    end
end


nMax=1e8; %maximum number of samples allowed.
out_param.alpha = alpha; %save the input parameters to a structure
out_param.inflate = inflate;
out_param.nSig = nSig;
tstart = tic; %start the clock 

%check to see input YYrand.
if isstruct(YYrand)
    if isfield(YYrand,'Yrand')
        Yrand=YYrand.Yrand;
    else
        Yrand = @(n) rand(n,1); %random number generator
    end
    if isfield(YYrand,'q')
        q=round(YYrand.q);
    else
        q=1;   
    end
else
    Yrand=YYrand;
    q=1;
end

r = size(Yrand(1),2);
p = max(0,r - q);
val = Yrand(nSig);

if p==0 && q==1
    Yval = val(:,1);
    samplevar = var(Yval);
    YY=Yval;
else 
        val = Yrand(nSig);
        meanVal=mean(val);
        A=bsxfun(@minus, val, meanVal);
        C=[ones(q,1); zeros(p,1)];
        [U,S,V]=svd([A; C'],0);
        Sdiag = diag(S);
        U2=U(end,:);
        y=U2'/(U2*U2');
        beta=V*(y./Sdiag);
        meanX=meanVal(:,q+1:end);
        meanX=[zeros(q,1); meanX'];
        YY = bsxfun(@minus, val, meanX')*beta;
        samplevar = var(YY);
end
    
out_param.sig0 = sqrt(samplevar); %standard deviation
sig0up = out_param.inflate.*out_param.sig0; %upper bound on the standard deviation
out_param.hmu0 = mean(YY);

nmu =  max(nSig,ceil((-gail.stdnorminv(alpha/2)*sig0up/max(absTol,relTol*abs(out_param.hmu0))).^2)); 
%nmu =  ceil((-gail.stdnorminv(alpha/2)*sig0up/max(absTol,relTol*abs(out_param.hmu0))).^2); 
   %number of samples needed for mean
if nmu > nMax %don't exceed sample budget
   warning(['The algorithm wants to use nmu = ' int2str(nmu) ...
      ', which is too big. Using ' int2str(nMax) ' instead.']) 
   nmu = nMax;
end

YY = Yrand(nmu);   
if p > 0 || q > 1
  YY(:,q+1:end) = bsxfun(@minus, YY(:,q+1:end), mean(YY(:,q+1:end),1));
  YY = YY*beta;
    %YY = val*beta;
end

hmu = mean(YY); %estimated mean
out_param.ntot = nSig+nmu; %total samples required
out_param.time = toc(tstart); %elapsed time
end


