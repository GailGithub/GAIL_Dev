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


%% Periodization transform
ptransform = 'C1sin'; %default option

if strcmp(ptransform,'Baker')
    f=@(x) f(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(ptransform,'C0')
    f=@(x) f(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(ptransform,'C1')
    f=@(x) f(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(ptransform,'C1sin')
    f=@(x) f(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
end




fx = f(x);
%% precompute the Bernoulli polynominal values to speedup the computation
out.aMLE(nn,1) = 0;
out_aMLE(nn,1) = 0;
out_costMLE(nn,1) = 0;
muhat(nn,1) = 0;
out_disc2(nn,1) = 0;
out_ErrBd(nn,1) = 0;

BernPolynOrder = 2;
%BernPolynX = zeros(size(x));

if strcmp(whKer,'Fourier')
    xdiff = abs(x - x(1,:));
    if BernPolynOrder==2
        BernPolynX = (6*xdiff.*(xdiff - 1) + 1)/6;
    elseif BernPolynOrder==4
        %BernPolynX = (30*(xdiff.^2.*(xdiff.^2 - 2*xdiff + 1)) - 1)/30;
        %BernPolynX = xdiff.^4 - 2*xdiff.^3 + xdiff.^2 -1/30;
        BernPolynX = (xdiff.^2).*(xdiff.^2 - 2*xdiff + 1) -(1.0/30);
    elseif BernPolynOrder==6
        BernPolynX = x.^6 - 3*x.^5 + (5/2)*x.^4 - (1/2)*x.^2 + (1/42);
    else
        fprintf('Error: Bernoulli polynomial order %d not implemented', BernPolynOrder);
        return;
    end
end

if  false %strcmp(whKer,'Fourier')
    %% plot MLEKernel cost function
    lnTheta = -3:0.05:0;
    costMLE = zeros(nn, numel(lnTheta));
    tic
    for ii = 1:nn
       nii = nvec(ii)
       xr = x(1:nii,:); fxr = fx(1:nii); BernPolynXr = BernPolynX(1:nii,:);
       %[~,w2] = sort(xr(:,1)); %w2 = bitrevorder(w2);
       %xr = xr(w2,:); fxr = fxr(w2); BernPolynXr = BernPolynXr(w2,:); %recorder to get symmetric matrix
       %par
       parfor k=1:numel(lnTheta)
           costMLE(ii,k) = MLEKernel(exp(lnTheta(k)),xr,fxr,whKer,domain,BernPolynXr,BernPolynOrder);
       end
    end
    toc
    figure; semilogx(exp(lnTheta),costMLE); lgd = legend(string(nvec),'location','north'); axis tight
    title(lgd,'Sample Size, \(n\)'); legend boxoff
    xlabel('Shape param, \(\theta\)')
    ylabel('MLE Cost, \( \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
    title(sprintf('MLE Cost function d=%d Bernoulli=%d PeriodTx=%s', d, BernPolynOrder, ptransform));
    [minVal,Index] = min(costMLE,[],2);
    hold on; plot(exp(lnTheta(Index)),minVal, '.');

    %print -depsc MLE_cost.eps
end

tstart = tic; %start the clock
parfor ii = 1:nn
   nii = nvec(ii);
   
   xr = x(1:nii,:); fxr = fx(1:nii); BernPolynXr = 0; % default Initial values
   if false % disable %strcmp(whKer,'Fourier')
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
   [lnaMLE, fval] = fminbnd(@(lna) ...
      MLEKernel(exp(lna),xr,fxr,whKer,domain,BernPolynXr,BernPolynOrder), ...
      ax,bx,optimset('TolX',1e-2));
   aMLE = exp(lnaMLE);

   %% Use the optimal \theta
   out_aMLE(ii) = aMLE
   out_costMLE(ii) = fval;
   [K,kvec,k0] = kernelFun(xr,whKer,aMLE,domain,BernPolynXr,BernPolynOrder);

   if strcmp(whKer,'Fourier')
       cn = K;

       fy = abs(fft_DIT(fxr));
       eigval = fft_DIT(cn');
       if any(eigval==0)
         fprintf('Zero eigval in estimating Kinvy \n')
       end
       eigvaln = eigval + (eigval==0)*eps;
       Kinvy = ifft(fy./eigvaln);

       %% compute the approximate mu
       muhat(ii) = sum(fxr)/sum(cn);
       muhat_old = real(kvec'*Kinvy);
       abs(muhat_old - muhat(ii))
   else
       Kinv = pinv(K);
       %w = Kinv*kvec;
       Kinvy = Kinv*fx(1:nii);
       
       %% compute the approximate mu
       muhat(ii) = real(kvec'*Kinvy);
   end

   %out_Kcond(ii) = cond(K);

   %% compute the discriminant
   if strcmp(whKer,'Fourier') == true
       %disc2_old = k0 - kvec'*ifft(fft_DIT(kvec)./fft_DIT(cn'));
       disc2 = k0 - nii/sum(cn);
       ffy = abs(fft_DIT(fxr));
       val2 = sum((ffy(eigval~=0).^2)./abs(eigval(eigval~=0))); %/nii;
       out_ErrBd(ii) = 2.58*sqrt(disc2*(val2))/nii;
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
       out_ErrBd(ii) = 2.58*sqrt(disc2*(fxr'*Kinvy)/nii);
   end
   out_disc2(ii) = disc2;
   
   
end

if strcmp(whKer,'Fourier') == true
    plot(out_aMLE,out_costMLE, 'c+');
end

out.aMLE = out_aMLE;
out.disc2 = out_disc2;
out.ErrBd = abs(out_ErrBd);
out.time = toc(tstart)
out.BernPolynOrder = BernPolynOrder;
out.ptransform = ptransform;

end


function y = fft_DIT( y )
nmmin = log2(length(y));
%y = bitrevorder(y);
for l=0:nmmin-1
    nl=2^l;
    nmminlm1=2^(nmmin-l-1);
    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
    coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
    coefv=repmat(coef,nmminlm1,1);
    evenval=y(ptind);
    oddval=y(~ptind);
    y(ptind)=(evenval+coefv.*oddval)/2;
    y(~ptind)=(evenval-coefv.*oddval)/2;
end
end
