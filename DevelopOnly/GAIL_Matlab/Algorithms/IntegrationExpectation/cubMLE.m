function [muhat,out]=cubMLE(f,nvec,domain,whSample,whKer,powerFuncMethod)
%CUBMLE Monte Carlo method to estimate the mean of a random variable
%
%   tmu = cubML(Yrand,absTol,relTol,alpha,nSig,inflate) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples may be of one of several kinds.  The default values are n=2^10 and
%   d = 1 Input f is a function handle that accepts an n x d matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%   powerFuncMethod: 
%      'cauchy' : using the Cauchy interlacing theorem
%      'thompson' : using the technique form R.C.Thompson paper
%   
% This is a heuristic algorithm based on a Central Limit Theorem
% approximation
if nargin < 7
    if ~exist('powerFuncMethod','var') || isempty(powerFuncMethod)
        powerFuncMethod='cauchy';  %technique to compute power function%
    end
    if nargin < 6
       if nargin < 5;
          whKer = 'Mat1'; %type of kernel
          if nargin < 4
             whSample = 'Sobol'; %type of sampling, scrambled Sobol
             if nargin < 3
                domain = [0;1]; %dimension
                if nargin < 2
                   nvec = 2^10; %number of samples
                   if nargin < 1
                      f = @(x) x.^2; %function
                   end
                end
             end
          end
       end
    end
end

tstart = tic; %start the clock
nmax = max(nvec);
nn = numel(nvec);
d = size(domain,2);
if strcmp(whSample,'Sobol')
   x = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
end
fx = f(x);
out.aMLE(nn,1) = 0;
muhat(nn,1) = 0;
for ii = 1:nn
   nii = nvec(ii);
   
   % Find the optimal shape parameter
   lnaMLE = fminbnd(@(lna) ...
      MLEKernel(exp(lna),x(1:nii,:),fx(1:nii),whKer,domain), ...
      -5,5,optimset('TolX',1e-2));
   aMLE = exp(lnaMLE);
   out.aMLE(ii) = aMLE;
   
   [K,kvec,k0] = kernelFun(x(1:nii,:),whKer,aMLE);
   Kinv = pinv(K);
   %w = Kinv*kvec;
   Kinvy = Kinv*fx(1:nii);
   
   % compute the approximate mu 
   muhat(ii) = kvec'*Kinvy;
   
   if strcmp(powerFuncMethod, 'cauchy')
       eigK = eig(K);
       eigKaug = eig([k0 kvec'; kvec K]);
       disc2 = exp(sum(log(eigKaug(1:nii)) - log(eigK)))*eigKaug(end);
   else
       % from R.C.Thompson paper
       [eigVecKaug, eigValKaug] = eig([k0 kvec'; kvec K]);
       if isvector(eigValKaug) == false
         eigValKaug = diag(eigValKaug);
       end
       uii = eigVecKaug(:,1);
       disc2 = 1/sum(uii.^2 ./ eigValKaug);
   end
   out.ErrBd(ii) = 2.58*sqrt(disc2*(fx(1:nii)'*Kinvy)/nii);
end  
out.time = toc(tstart);
