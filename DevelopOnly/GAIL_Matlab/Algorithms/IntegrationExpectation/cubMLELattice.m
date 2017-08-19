function [muhat,out]=cubMLELattice(f,d,absTol,relTol,order,ptransform, ...
  testAll,figSavePath,fName,arbMean)
%CUBMLE Monte Carlo method to estimate the mean of a random variable
%
%   tmu = cubMLELattice(f,absTol,relTol,alpha,nSig,inflate) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples may be of one of several kinds.  The default values are n=2^10 and
%   d = 1 Input f is a function handle that accepts an n x d matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%
%
% This is a heuristic algorithm based on a Central Limit Theorem
% approximation

if nargin < 6
  ptransform = 'Baker';
  if nargin < 5
    order = 2; %type of sampling, scrambled Sobol
    if nargin < 4
      relTol = 0;
      if nargin < 3
        absTol = 0.01;
        if nargin < 2
          d = 1; %dimension
          if nargin < 1
            f = @(x) x.^2; %function
          end
        end
      end
    end
  end
end

if ~exist('testAll','var')
  testAll = false;
end

debugEnable=false;

% comment this line of code to use GPU for computations
gpuArray = @(x) x;   gather = @(x) x;

tstart = tic; %start the clock

% define min and max number of points allowed in the automatic cubature
mmin = 10;
mmax = 23;
mvec = mmin:mmax;
numM = length(mvec);

% variables to save debug info
errorBdAll = zeros(length(mvec),1);
muhatAll = zeros(length(mvec),1);
aMLEAll = zeros(length(mvec),1);
timeAll = zeros(length(mvec),1);
dscAll = zeros(length(mvec),1);
s_All = zeros(length(mvec),1);

% pick a random value to apply as shift
shift = rand(1,d);

% apply periodization transformation to the function
ff = doPeriodTx(f, ptransform);


%% plot MLE loss function
if exist('fName','var') && exist('figSavePath','var')
  plotMLE_Loss(ff, mvec, figSavePath, fName, d, order, ptransform,arbMean,visiblePlot)
end


%% Iteratively find the number of points required for the cubature to meet
% the error threshold
for ii = 1:numM
  tstart_iter = tic; 
  m = mvec(ii);
  n = 2^m;
  
  %Update function values
  %% Efficient FFT computation algorithm
  if ii == 1
    % in the first iteration compute the full FFT
    xun = simple_lattice_gen(n,d,true);
    x = mod(bsxfun(@plus,xun,shift),1);  % shifted
    
    ftildeNew=gpuArray(ff(x)); %evaluate integrand
    
    % Compute initial FFT
    for l=0:mmin-1
      nl=2^l;
      nmminlm1=2^(mmin-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
      coefv=repmat(coef,nmminlm1,1);
      evenval=ftildeNew(ptind);
      oddval=ftildeNew(~ptind);
      ftildeNew(ptind)=(evenval+coefv.*oddval);
      ftildeNew(~ptind)=(evenval-coefv.*oddval);
    end
    
  else
    xunnew = simple_lattice_gen(n,d,false);
    xnew = mod(bsxfun(@plus,xunnew,shift),1);
    xun = [xun;xunnew];
    x = [x;xnew];
    
    mnext=m-1;
    
    ftildeNextNew=gpuArray(ff(xnew));  % initialize for inplace computation
    
    if debugEnable==true
      if any(isnan(ftildeNextNew)) || any(isinf(ftildeNextNew))
        fprintf('ftildeNextNew NaN \n');
      end
    end
    
    % Compute initial FFT on next set of new points
    for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
      coefv=repmat(coef,nmminlm1,1);
      evenval=ftildeNextNew(ptind);
      oddval=ftildeNextNew(~ptind);
      ftildeNextNew(ptind)=(evenval+coefv.*oddval);
      ftildeNextNew(~ptind)=(evenval-coefv.*oddval);
      
      if debugEnable==true
        if any(isnan(ftildeNextNew)) || any(isinf(ftildeNextNew))
          fprintf('ftildeNextNew NaN \n');
        end
      end
      
    end
    
    % combine the previous batch and new batch to get FFT on all points
    ftildeNew=[ftildeNew;ftildeNextNew];
    nl=2^mnext;
    ptind=[true(nl,1); false(nl,1)];
    coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
    coefv=repmat(coef,nmminlm1,1);
    evenval=ftildeNew(ptind);
    oddval=ftildeNew(~ptind);
    ftildeNew(ptind)=(evenval+coefv.*oddval);
    ftildeNew(~ptind)=(evenval-coefv.*oddval);
    
  end
  
  ftilde = ftildeNew;
  % ftilde(1) = sum(ff(x)); % correction to avoid round off error
  br_xun = bitrevorder(gpuArray(xun));
  
  %Compute MLE parameter
  lnaMLE = fminbnd(@(lna) ...
    MLEKernel(exp(lna),br_xun,ftilde,order,arbMean), ...
    -5,0,optimset('TolX',1e-2)); % -5,5
  aMLE = exp(lnaMLE);
  [loss,Ktilde,KtildeSq,RKHSnorm,K] = MLEKernel(aMLE,br_xun,ftilde,order,arbMean);
  wt = 1/Ktilde(1);
  
  %Check error criterion
  DSC = abs(1/n - wt);
  
  % store the debug information
  dscAll(ii) = DSC;
  s_All(ii) = RKHSnorm;
  
  out.ErrBd = 2.58*sqrt(DSC*RKHSnorm);
  if arbMean==true % zero mean case
    muhat = ftilde(1)/n;
  else % non zero mean case
    muhat = ftilde(1)/Ktilde(1);
  end
  muminus = muhat - out.ErrBd;
  muplus = muhat + out.ErrBd;
  muhatAll(ii) = muhat;
  errorBdAll(ii) = out.ErrBd;
  aMLEAll(ii) = aMLE;
  
  if 2*out.ErrBd <= ...
      max(absTol,relTol*abs(muminus)) + max(absTol,relTol*abs(muplus))
    
    if errorBdAll(ii)==0
      errorBdAll(ii) = eps;
    end
    
    % if testAll flag is set, run for for all 'n' values to compute error
    % used for error plotting
    if testAll==false
      break
    end
  end
  
  timeAll(ii) = toc(tstart_iter);
end

out.n = n;
out.time = toc(tstart);
out.ErrBdAll = errorBdAll;
out.muhatAll = muhatAll;
out.mvec = mvec;
out.aMLEAll = aMLEAll;
out.timeAll = timeAll;
out.s_All = s_All;
out.dscAll = dscAll;

% convert from gpu memory to local
muhat=gather(muhat);
out=gather(out);

end

function [K, Ktilde] = kernel(xun,order,a)

constMult = -(-1)^(order/2)*(2*pi)^order/factorial(order);
if order == 2
  bernPloy = @(x)(-x.*(1-x) + 1/6);
elseif order == 4
  %bernPloy = @(x)((x.^2).*(x.*(x-2) +1) - 1/30);
  bernPloy = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
else
  error('Bernoulli order not implemented !');
end
K = prod(1 + a*constMult*bernPloy(xun),2);

% matlab's builtin fft is much faster and accurate
Ktilde = real(fft(K));

if sum(K)==length(K) || Ktilde(1)==length(K)
  %fprintf('debug');
end
end

function [loss,Ktilde,KtildeSq,RKHSnorm,K] = MLEKernel(a,xun,ftilde,...
  order,useArbMean)

roundOffErrOptimize = false;
n = length(ftilde);
if order==4
  if n>(2^18) && roundOffErrOptimize
    [K, Ktilde] = kernel(xun,order,a);
    [K2, Ktilde2] = kernel(xun,order/2,sqrt(a));
    
    KtildeSq = (Ktilde2/sqrt(n));
  else
    [K, Ktilde] = kernel(xun,order,a);
    KtildeSq = sqrt(Ktilde);
  end
elseif order==2
  [K, Ktilde] = kernel(xun,order,a);
  KtildeSq = sqrt(Ktilde);
else
  error('Unsupported Bernoulli polyn order !');
end

Ktilde = abs(Ktilde);  % remove any negative values
ftilde = abs(ftilde);

if any(KtildeSq==0)
  % fprintf('Ktilde has zero vals \n');
end

%RKHSnorm = mean(abs(ftilde).^2./Ktilde);

% temp = (abs(ftilde(KtildeSq~=0))./(KtildeSq(KtildeSq~=0))).^2 ;
temp = (abs(ftilde(Ktilde~=0).^2)./(Ktilde(Ktilde~=0))) ;

if useArbMean==true
  RKHSnorm = sum(temp(2:end))/numel(Ktilde);
  temp_1 = sum(temp(2:end));
else
  RKHSnorm = sum(temp)/numel(Ktilde);
  temp_1 = sum(temp);
end
RKHSnormSq = sqrt(RKHSnorm);

if isnan(RKHSnormSq)
  fprintf('RKHSnormSq NaN \n');
end
% loss = mean(2*log(KtildeSq)) + log(RKHSnorm);
% loss = mean(log(Ktilde)) + log(RKHSnorm);
loss = sum(log(Ktilde)) + numel(Ktilde)*log(temp_1);

if ~isreal(Ktilde)
  fprintf('Ktilde has complex vals \n');
end

if ~isreal(temp_1)
  fprintf('temp_1 has complex vals \n');
end

if isnan(loss)
  fprintf('loss NaN \n');
end

a;
%ftilde(1)/Ktilde(1);
end

function f = doPeriodTx(fInput, ptransform)

if strcmp(ptransform,'Baker')
  f=@(x) fInput(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(ptransform,'C0')
  f=@(x) fInput(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(ptransform,'C1')
  f=@(x) fInput(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(ptransform,'C1sin')
  f=@(x) fInput(x-sin(2*pi*x)/(2*pi)).*prod(2*sin(pi*x).^2,2); % Sidi C^1 transform
elseif strcmp(ptransform,'C2sin')
  psi3 = @(t) (8-9*cos(pi*t)+cos(3*pi*t))/16; psi3_1 = @(t) (9*sin(pi*t)*pi- sin(3*pi*t)*3*pi)/16;
  f=@(x) fInput(psi3(x)).*prod(psi3_1(x),2);
elseif strcmp(ptransform,'C3sin')
  psi4 = @(t) (12*pi*t-8*sin(2*pi*t)+sin(4*pi*t))/(12*pi); psi4_1 = @(t) (12*pi-8*cos(2*pi*t)*2*pi+sin(4*pi*t)*4*pi)/(12*pi);
  f=@(x) fInput(psi4(x)).*prod(psi4_1(x),2);
elseif strcmp(ptransform,'none')
  % do nothing
  f=@(x) fInput(x);
else
  error('Error: Periodization transform %s not implemented', ptransform);
end

end

function xlat = simple_lattice_gen(n,d,firstBatch)
if d<=10
  % this gives best accuracy
  z = [1, 364981, 245389, 97823, 488939, 62609, 400749, 385317, 21281, 223487]; % generator from Hickernell's paper
  %z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599]; %generator
else
  z = [1 182667 302247 433461 160317 94461 481331 252345 358305 221771 48157 489023 438503 399693 200585 169833 308325 247437 281713 424209 244841 205461 336811 359375 86263 370621 422443 284811 231547 360239 505287 355195 52937 344561 286935 312429 513879 171905 50603 441451 164379 139609 371213 152351 138607 441127 157037 510073 281681 380297 208143 497641 482925 233389 238553 121499 137783 463115 168681 70699];
end

z = z(1:d);

if false
  % this is very slow
  if firstBatch==true
    brIndices = bitrevorder((0:1/n:1-1/n));
  else
    brIndices = bitrevorder((1/n:2/n:1-1/n));
  end
else
  nmax = n;
  nmin = 1 + n/2;
  if firstBatch==true
    nmin = 1;
  end
  nelem=nmax-nmin+1;
  
  if firstBatch==true
    brIndices=vdc(nelem)';
  else
    brIndices=vdc(nelem)'+1/(2*(nmin-1));
  end
end
xlat = mod(bsxfun(@times,brIndices',z),1);  % unshifted

end

% Van der Corput sequence in base 2
function q = vdc(n)
if n>1
  k=log2(n); % We compute the VDC seq part by part of power 2 size
  q=zeros(2^k,1);
  for l=0:k-1
    nl=2^l;
    kk=2^(k-l-1);
    ptind=repmat([false(nl,1);true(nl,1)],kk,1);
    q(ptind)=q(ptind)+1/2^(l+1);
  end
else
  q=0;
end
end


% fft with deimation in time i.e, input is already in 'bitrevorder'
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


function minTheta = plotMLE_Loss(ff, mvec, fullPath, fName, d, order, ptransform,useArbMean,visiblePlot)

n = 2.^mvec(end);
xun = simple_lattice_gen(n,d,true);
fx = ff(xun);  % Note: periodization transform already applied
numM = length(mvec);

%% plot MLEKernel cost function
lnTheta = -3:0.1:5;
% fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
plotFileName = sprintf('%s%s Cost d_%d bernoulli_%d Period_%s.eps', fullPath, fName, d, order, ptransform)

costMLE = zeros(numM,numel(lnTheta));
tic

for ii = 1:numM
  nii = 2^mvec(ii)
  
  eigvalK = zeros(numel(lnTheta),nii);
  ftilde = fft(bitrevorder(fx(1:nii))); %/nii;
  
  tic
  %par
  parfor k=1:numel(lnTheta)
    [costMLE(ii,k),eigvalK(k,:)] = MLEKernel(exp(lnTheta(k)),xun(1:nii,:),ftilde,order,useArbMean);
  end
  toc
end

toc

if exist('visiblePlot','var') && visiblePlot==false
  hFigCost = figure('visible','off');
else
  hFigCost = figure();
end
semilogx(exp(lnTheta),real(costMLE)); %lgd = legend(string(nvec),'location','north'); axis tight
set(hFigCost, 'units', 'inches', 'Position', [4 4 10 7])
%title(lgd,'Sample Size, \(n\)'); legend boxoff
xlabel('Shape param, \(\theta\)')
ylabel('MLE Cost, \( \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
axis tight;
if useArbMean
  mType = '\(m \neq 0\)';
else
  mType = '\(m = 0\)';
end
title(sprintf('%s d=%d r=%d Tx=%s %s', fName, d, order, ptransform, mType));
[minVal,Index] = min(real(costMLE),[],2);
hold on; plot(exp(lnTheta(Index)),minVal, '.');
saveas(hFigCost, plotFileName)

minTheta = exp(lnTheta(Index));

end
