function [muhat,out]=cubBayesMLE_Matern_g(varargin)
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

% f,nvec,domain,whSample,whKer,powerFuncMethod
default_args.f = @(x) x.^2; %function
default_args.nvec = 2^10; %number of samples
default_args.whSample = 'Sobol'; %type of sampling, scrambled Sobol
default_args.whKer = 'Mat1'; %type of kernel
default_args.powerFuncMethod = 'cauchy';
s_args = parse_args(default_args, varargin{:});

f = s_args.f;
whKer = s_args.whKer;
domain = repmat([0;1], 1,s_args.dim);

mvec = 5:1:13;  %15; % 2^14 comp takes time > 5 hours
nvec = 2.^mvec;
tstart = tic; %start the clock
nmax = max(nvec);
nn = numel(nvec);
d = size(domain,2);
if strcmp(s_args.whSample,'Sobol')
  x = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
end
fx = f(x);
out.aMLE(nn,1) = 0;
muhatVec(nn,1) = 0;
muhat = 0;
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
  %Kinvy = Kinv*fx(1:nii);
  Kinvy = K\fx(1:nii);
  
%   if sum(abs(Kinvy_-Kinvy)) > 1
%     fprintf('Error too big');
%   end
  
  % compute the approximate mu
  muhatVec(ii) = kvec'*Kinvy;
  
  if strcmp(s_args.powerFuncMethod, 'cauchy')
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
  out.ErrBdVec(ii) = 2.58*sqrt(disc2*(fx(1:nii)'*Kinvy)/nii);
  muminus = muhatVec(ii) - out.ErrBdVec(ii);
  muplus = muhatVec(ii) + out.ErrBdVec(ii);
  
  if 2*out.ErrBdVec(ii) <= ...
      max(s_args.absTol,s_args.relTol*abs(muminus)) + max(s_args.absTol,s_args.relTol*abs(muplus))
    muhat = muhatVec(ii);
    break
  end
end
out.time = toc(tstart);
out.n = nii;
out.absTol = s_args.absTol;
out.reltol = s_args.relTol;
end

function val = MLEKernel(shape,x,y,whKer,domain)
% The function MLEKERNEL is the loss function to be minimized to obtain the
% shape parameter |shape| that defines the kernel.  It uses data |(x,y)|
% and the function |kernelFun|.
if nargin < 4
  whKer = 'Mat1';
end
nx = size(x,1);
K = kernelFun(x,whKer,shape,domain);
[eigvec,eigval] = eig(K,'vector');
Vty = eigvec'*y;
%val = sum(log(eigval))/nx + Vty'*(Vty./eigval);
val = sum(log(eigval))/nx + log(Vty'*(Vty./eigval));  %corrected
end

function [K,kvec,k0] = kernelFun(x,whKer,shape,domain)
[nx,d] = size(x);
if nargin < 4
  domain = [zeros(1,d); ones(1,d)];
  if nargin < 3
    shape = 1;
    if nargin < 2
      whKer = 'sqExp';
    end
  end
end
K = ones(nx);
if strcmp(whKer,'sqExp') %al  l wrong
  kvec = ones(nx,1)*(sqrt(pi)/(2*shape))^d;
  for k = 1:d
    K = K.*exp(-(shape*bsxfun(@minus,x(:,k),x(:,k)')).^2);
    kvec = kvec.*(erf(shape*x(:,k)) + erf(shape*(1 - x(:,k))));
  end
elseif strcmp(whKer,'Mat1')
  diffdom = diff(domain,1,1);
  shdiffdom = shape*diffdom;
  k0 = prod((- 6 + 4*shdiffdom +exp(-shdiffdom).*(6 + 2*shdiffdom))./shdiffdom.^2);
  kvec = ones(nx,1)*(2^d/prod(shdiffdom));
  for k = 1:d
    tempa = shape*abs(bsxfun(@minus,x(:,k),x(:,k)'));
    K = K.*exp(-tempa).*(1 + tempa);
    tempb = shape*(x(:,k)-domain(1,k));
    tempc = shape*(domain(2,k) - x(:,k));
    kvec = kvec.*(2 - exp(-tempc).*(1+tempc/2) ...
      - exp(-tempb).*(1+tempb/2));
  end
end

end

function s_args = parse_args(s_args, varargin)

iStart = 1;

% parse each input argument passed
if nargin >= iStart
  wh = find(strcmp(varargin(iStart:end),'f'));
  if ~isempty(wh), s_args.f = varargin{wh+iStart}; end
  wh = find(strcmp(varargin(iStart:end),'dim'));
  if ~isempty(wh), s_args.dim = varargin{wh+iStart}; end
  wh = find(strcmp(varargin(iStart:end),'absTol'));
  if ~isempty(wh), s_args.absTol = varargin{wh+iStart}; end
  wh = find(strcmp(varargin(iStart:end),'relTol'));
  if ~isempty(wh), s_args.relTol = varargin{wh+iStart}; end
  wh = find(strcmp(varargin(iStart:end),'stopAtTol'));
  if ~isempty(wh), s_args.stopAtTol = varargin{wh+iStart}; end
  wh = find(strcmp(varargin(iStart:end),'fName'));
  if ~isempty(wh), s_args.fName = varargin{wh+iStart}; end
end

end
