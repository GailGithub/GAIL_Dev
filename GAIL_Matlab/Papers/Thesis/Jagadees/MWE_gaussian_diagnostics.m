%
% Minimum working example to test Gaussian diagnostics idea
%
function MWE_gaussian_diagnostics()


format short
close all
fNames = {'ExpCos','MVN','Keister','rand'};
ptransforms = {'C1','C1sin', 'none'};
fName = fNames{2};
ptransform = ptransforms{2}; 
whKer = 'Bern';  %'Exp'; %

npts = 2^10;  % max 14
dim = 1;
shift = rand(1,dim);

[~,xlat_] = cubBayesLattice_g.simple_lattice_gen(npts,dim,shift,true);

if strcmp(fName,'ExpCos')
  integrand = @(x) exp(sum(cos(2*pi*x), 2));
elseif strcmp(fName,'MVN')
  dim=2; absTol=1e-3; relTol=1e-2; 
  C = [4 1 1; 0 1 0.5; 0 0 0.25]; MVNParams.Cov = C'*C; MVNParams.C = C;
  MVNParams.a = [-6 -2 -2]; MVNParams.b = [5 2 1]; MVNParams.mu = 0;
  MVNParams.CovProp.C = chol(MVNParams.Cov)';
  exactF = 0.676337324357787;
  integrand = @(t) GenzFunc(t,MVNParams);
elseif strcmp(fName, 'Keister')
  integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8
else
  integrand = @(x) f_rand(x, 1);
end

rng(202326) % initialize random number generator for reproducability

N = 2^(15);
f_c = randn(N, 1);
f_s = randn(N, 1);
f_0 = randn(1, 1);
  
% gaussian random function
function fval = f_rand(xpts,theta)

  kvec = (1:N)';
  argx = @(x) 2*pi*bsxfun(@times, kvec, x);
  f_c_ = @(x)(f_c./kvec).*cos(argx(x));
  f_s_ = @(x)(f_s./kvec).*sin(argx(x));
  f_ran = @(x,theta) prod((f_0 + theta * sum(f_c_(x) + f_s_(x) )), 2) ;
  [n,~]=size(xpts);
  fval=zeros(n,1);
  for i=1:n
    fval(i) = f_ran(xpts(i,:),theta);
  end
end

% bernPoly = @(x)(-x.*(1-x) + 1/6);

% fval = bernPoly(0.5);
% fval_ = bernoulli_series(0.5,2);

% fval = bernPoly(xlat_);
% fval_ = cubBayesLattice_g.bernoulli_series(xlat_,2);
% sum(sum(abs(fval - fval_)))

ff = f_rand(xlat_, 1);
[H_,pValue_,W_] = swtest(ff)
if H_==1, Hval='false'; else, Hval='true'; end
fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3e, W=%1.3f, n=%d\n', ...
  Hval, pValue_, W_, npts);
      
% integrand = @(x) prod(bernPoly(x), 2);

inputArgs = {'f',integrand,'fName',fName,'dim',dim, 'order',1, ....
  'ptransform',ptransform, 'stopAtTol',true, 'stopCriterion','MLE', ...
  'arbMean',true, 'visiblePlot',true};

obj=cubBayesLattice_g(inputArgs{:});

[minTheta, figH] = obj.plotObjectiveFunc();
% build filename with path to store the plot
plotFileName = sprintf('%s%s_Cost_d%d_r%d_%s_%s.png',...
  obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform, obj.stopCriterion);
saveas(figH, plotFileName)




% integrand = @(x) x.^2;
integrand_p = cubBayesLattice_g.doPeriodTx(integrand, ptransform);

y = integrand_p(xlat_);
ftilde = fft(y);
ftilde(1) = 0;  % ((f - m_\MLE I)
if dim==1
  figure; scatter(xlat_, y, 10)
  title('Integrand')
end


obj.f = integrand_p;
if strcmp(whKer, 'Bern')
  obj.kernType = 1;
else
  obj.kernType = 2;
end
obj.dim = dim;
obj.ptransform = ptransform;
obj.avoidCancelError = true;
obj.debugEnable = false;
obj.order = 1;
obj.arbMean = true;
obj.visiblePlot = true;
obj.fName = fName;
obj.figSavePath = '';

%minTheta = plotObjectiveFunc(obj)

% [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj,0.5,xlat_,ftilde);
lnaMLE = fminbnd(@(lna) ...
  ObjectiveFunction(obj, exp(lna),xlat_,(ftilde)), ...
  -5,5,optimset('TolX',1e-2));
aMLE = exp(lnaMLE);

lnaMLE_opt = fminsearch(@(lna) ...
  ObjectiveFunction(obj, exp(lna),xlat_,(ftilde)), ...
  0,optimset('TolX',1e-2));
thetaOpt = exp(lnaMLE_opt);


% bVec = [0.3]; %[0.3:0.1:0.9];  %
bVec = [obj.order];
%for sh=[0.1 0.5 0.9 5]  %0.001 0.01
for b=bVec
  % lambda = kernel(xlat_,whKer,thetaOpt,b,obj.avoidCancelError);
  debugEnable = false;
  [lambda] = cubBayesLattice_g.kernel(xlat_,obj.order,thetaOpt,obj.avoidCancelError,...
          obj.kernType,debugEnable);
      
  w_ftilde = real(ifft(ftilde./sqrt(abs(real(lambda)))));
  % figure; plot(abs(real(lambda)), '.')
  % figure; plot(abs(real(ftilde)), '.')
  % figure; plot(w_ftilde, '.')
  hFigNormplot = figure();
  set(hFigNormplot,'defaultaxesfontsize',16, ...
    'defaulttextfontsize',16, ... %make font larger
    'defaultLineLineWidth',0.75, 'defaultLineMarkerSize',8)
  normplot(w_ftilde)
  
  % Shapiro-Wilk test
  [H, pValue, W] = swtest(w_ftilde);
  Hval='true';
  if H==true
    Hval='false';
  end
  fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3f, W=%1.3f\n', Hval, pValue, W);
  
  title(sprintf('%s n=%d Tx=%s b=%1.2f theta=%1.2f', ...
    fName, npts, ptransform, b, thetaOpt))
end
%end

end

% function [Lambda, Lambda_ring] = kernel(xun,order,a,avoidCancelError,...
%         kernType,debugEnable)
%{
function [Lambda, Lambda_ring] = kernel(xun,whKer,a,b,avoidCancelError)

if strcmp(whKer,'Bern')
  b_order = b;
  constMult = -(-1)^(b_order/2)*((2*pi)^b_order)/factorial(b_order);
  if b_order==2
    kernelFunc = @(x)(-x.*(1-x) + 1/6);
  elseif b_order==4
    kernelFunc = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
  else
    error('Bernoulli order not implemented !');
  end
  
  
else
  %b = 0.7;
  kernelFunc = @(x) 2*b*((cos(2*pi*x)-b))./(1 + b^2 - 2*b*cos(2*pi*x));
  constMult = 1;
end

if avoidCancelError
  [C1m1, C1_alt] = cubBayesLattice_g.kernel_t(a*constMult, kernelFunc(xun));
  % eigenvalues must be real : Symmetric pos definite Kernel
  Lambda_ring = real(fft(C1m1));
  
  Lambda = Lambda_ring;
  Lambda(1) = Lambda_ring(1) + length(Lambda_ring);
else
  C1 = prod(1 + (a)*constMult*kernelFunc(xun),2);  %
  Lambda = real(fft(C1));
  Lambda_ring = nan;
end

end
%}

function [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj,a,xun,ftilde)

n = length(ftilde);
% [Lambda, Lambda_ring] = kernel(xun,obj.kernType,a,obj.order,...
%   obj.avoidCancelError);

debugEnable = false;
[Lambda, Lambda_ring] = cubBayesLattice_g.kernel(xun,obj.order,a,obj.avoidCancelError,...
        obj.kernType,debugEnable);
      
% compute RKHSnorm
temp = abs(ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0)) ;

% compute loss: MLE
if obj.arbMean==true
  RKHSnorm = sum(temp(2:end))/n;
  temp_1 = sum(temp(2:end));
else
  RKHSnorm = sum(temp)/n;
  temp_1 = sum(temp);
end

% ignore all zero eigenvalues
loss1 = sum(log(Lambda(Lambda~=0)))/n;
loss2 = log(temp_1);
loss = (loss1 + loss2);
fprintf('loss1 %1.3f loss2 %1.3f loss %1.3f a %1.3e\n', loss1, loss2, loss, a)

if obj.debugEnable
  cubMLESobol.alertMsg(RKHSnorm, 'Imag');
  cubMLESobol.alertMsg(loss1, 'Inf');
  cubMLESobol.alertMsg(loss2, 'Inf');
  cubMLESobol.alertMsg(loss, 'Inf', 'Imag', 'Nan');
  cubMLESobol.alertMsg(Lambda, 'Imag');
end
end


% plots the objective function for the MLE of theta
function minTheta = plotObjectiveFunc(obj)

obj.mmax = 20;
obj.mvec = 6:1:10;
n = 2^obj.mmax;
shift = rand(1,obj.dim);

[~,xun_] = cubBayesLattice_g.simple_lattice_gen(n,obj.dim,shift,true);

xpts = xun_;
fx = obj.f(xpts);  % No periodization transform required
numM = length(obj.mvec);

%% plot MLEKernel cost function
lnTheta = -10:0.2:10;

%fullPath = strcat(obj.figSavePath,'/',obj.fName,'/',obj.ptransform,'/');
plotFileName = sprintf('%s%s Cost d_%d order_%d Period_%s.png', ...
  obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform);
plotFileName  % just to display it

costMLE = zeros(numM,numel(lnTheta));
tstart=tic;

for iter = 1:numM
  nii = 2^obj.mvec(iter);
  nii
  
  eigvalK = zeros(numel(lnTheta),nii);
  %br_xpts = bitrevorder(xpts(1:nii, :));
  
  ftilde = fft(fx(1:nii));
  br_xpts = xpts(1:nii,:);
  
  tic
  %par
  for k=1:numel(lnTheta)
    [costMLE(iter,k),eigvalK(k,:)] = ObjectiveFunction(obj,exp(lnTheta(k)),...
      br_xpts,ftilde);
  end
  toc
end

toc(tstart)

if obj.visiblePlot==false
  hFigCost = figure('visible','off');
else
  hFigCost = figure();
end

if ~isreal(costMLE)
  fprintf('costMLE has complex values !! \n')
end

% semilogx : loglog
semilogx(exp(lnTheta),real(costMLE));
set(hFigCost, 'units', 'inches', 'Position', [0 0 13.5 11.5])
xlabel('Shape param, \(\theta\)')
ylabel('MLE Cost, \( \log \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
% ylabel('Log MLE Obj. fun.')
axis tight;
if obj.arbMean
  mType = '\(m \neq 0\)';  % arb mean
else
  mType = '\(m = 0\)';  % zero mean
end
%title(lgd,'Sample Size, \(n\)'); legend boxoff
title(sprintf('%s d=%d r=%d %s', obj.fName, obj.dim, ...
  obj.order, mType));
[minVal,Index] = min(real(costMLE),[],2);

% mark the min theta values found using fminbnd
minTheta = exp(lnTheta(Index));
hold on;
semilogx(minTheta,minVal, '.');
if exist('aMLEAll', 'var')
  semilogx(obj.aMLEAll,obj.lossMLEAll, '+');
end
temp = string(obj.mvec); temp=strcat('\(2^{',temp,'}\)');
temp(end+1) = '\(\theta_{min_{true}}\)';
if exist('aMLEAll', 'var')
  temp(end+1) = '\(\theta_{min_{est}}\)';
end

legend(temp,'location','best'); axis tight
saveas(hFigCost, plotFileName)
end



