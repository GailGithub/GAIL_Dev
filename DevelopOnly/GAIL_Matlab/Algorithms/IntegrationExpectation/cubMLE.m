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
        powerFuncMethod='Cauchy';  %technique to compute power function%
    end
    if nargin < 6
       if nargin < 5
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

nmax = max(nvec);
nn = numel(nvec);
d = size(domain,2);
if strcmp(whSample,'Sobol')
   x = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
end

%% Generate Rank-1 Lattice points
if strcmp(whSample,'Lattice1') 
    x = gail.lattice_gen(1,nmax,d);
end

fx = f(x);
%% precompute the Bernoulli polynominal values to speedup the computation
out.aMLE(nn,1) = 0;
out_aMLE(nn,1) = 0;
muhat(nn,1) = 0;
out_disc2(nn,1) = 0;
out_ErrBd(nn,1) = 0;
BernPolynOrder = 2;
BernPolynX = zeros(size(x));
if strcmp(whKer,'Fourier')
    xdiff = abs(x - x(1,:));
    if BernPolynOrder==2
        BernPolynX = (6*xdiff.*(xdiff - 1) + 1)/6;
    elseif BernPolynOrder==4
        BernPolynX = (30*(xdiff.^2.*(xdiff.^2 - 2*xdiff + 1)) - 1)/30;
    else
        fprintf('Bernoulli polynomial order not implemented');
    end
end

if  false %strcmp(whKer,'Fourier')
    %% plot MLEKernel cost function
    costX = -3:0.05:0;
    costMLE = zeros(nn, numel(costX));
    tic
    for ii = 1:nn
       nii = nvec(ii)
       xr = x(1:nii,:); fxr = fx(1:nii); BernPolynXr = BernPolynX(1:nii,:);
       [~,w2] = sort(xr(:,1));xr = xr(w2,:); fxr = fxr(w2); %recorder to get symmetric matrix
       BernPolynXr = BernPolynXr(w2,:);
       parfor k=1:numel(costX) 
           costMLE(ii,k) = MLEKernel(exp(costX(k)),xr,fxr,whKer,domain,BernPolynXr,BernPolynOrder);
       end
    end
    toc
    figure; semilogx(exp(costX),costMLE); lgd = legend(string(nvec),'location','north'); axis tight
    title(lgd,'Sample Size, \(n\)'); legend boxoff
    xlabel('Shape param, \(\theta\)')
    ylabel('MLE Cost, \( \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
    title('MLE Cost function');
    print -depsc MLE_cost.eps
end

tstart = tic; %start the clock
for ii = 1:nn
   nii = nvec(ii);
   
   xr = x(1:nii,:); fxr = fx(1:nii); BernPolynXr = 0; % default Initial values
   if strcmp(whKer,'Fourier')
       %% Reorder the x, fx to get symmetric circulant (Toeplitz) matrix
       [~,w2]=sort(xr(:,2)); xr = xr(w2,:);
       fxr = fxr(w2);
       BernPolynXr = BernPolynX(w2,:);
   end
   if strcmp(whKer,'Fourier')
       ax = -3; bx = 0; % limit theta within (0,1]
   else
       ax = -5; bx = 5;
   end
   
   %% Estimate optimal \theta
   lnaMLE = fminbnd(@(lna) ...
      MLEKernel(exp(lna),xr,fxr,whKer,domain,BernPolynXr,BernPolynOrder), ...
      ax,bx,optimset('TolX',1e-2));
   aMLE = exp(lnaMLE);
   
   %% Use the optimal \theta
   out_aMLE(ii) = aMLE
   [K,kvec,k0] = kernelFun(xr,whKer,aMLE,domain,BernPolynXr,BernPolynOrder);
   
   if strcmp(whKer,'Fourier')
       y = fxr; cn = K;
       Kinvy = ifft(fft(y)./fft(cn'));
   else
       Kinv = pinv(K);
       %w = Kinv*kvec;
       Kinvy = Kinv*fx(1:nii);
   end
   
   %% compute the approximate mu 
   muhat(ii) = kvec'*Kinvy
   %out_Kcond(ii) = cond(K);

   %% compute the discriminant
   if strcmp(whKer,'Fourier') == true
       disc2 = k0 - kvec'*ifft(fft(kvec)./fft(cn'));
   else
       if strcmp(powerFuncMethod, 'Cauchy')
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
           disc2 = 1.0/sum(uii.^2 ./ eigValKaug);
       end
   end
   out_disc2(ii) = disc2;
   out_ErrBd(ii) = 2.58*sqrt(disc2*(fxr'*Kinvy)/nii);
end  

out.aMLE = out_aMLE;
out.disc2 = out_disc2;
out.ErrBd = out_ErrBd;
out.time = toc(tstart)

end

