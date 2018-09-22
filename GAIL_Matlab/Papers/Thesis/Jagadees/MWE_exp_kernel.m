%
% Minimum working example to test the new exponentially converging kernel
%
function MWE_exp_kernel()
format short

npts = 2^12;  % max 14
dim = 2;
shift = rand(1,dim);
z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
            151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
            85729, 14597, 94813, 422013, 484367]; %generator      
z = z(1:dim);

bernPoly = @(x)(-x.*(1-x) + 1/6); kern = @(x) prod((1 + 0.9*bernPoly(x)), 2);
xlat_ = mod(bsxfun(@times,(0:1/npts:1-1/npts)',z),1); % unshifted in direct order
xlat = mod(bsxfun(@plus,xlat_,shift),1);  % shifted in direct order

if 0
  C1 = kern(xlat_);
  K = ones(npts);
    % tempa = abs(bsxfun(@minus,x(:,k),x(:,k)') );  % Lattice
  tempa = mod(bsxfun(@minus,xlat_,xlat_'), 1);  % Lattice
  kern = @(x) ((1 + 0.9*bernPoly(x)));
  K = kern(tempa);
  xx = (1:npts) - 1;
  indx = mod((xx'*xx), npts);
  V = exp(-sqrt(-1)*2*pi*(indx)/npts);
  eigK = fft(K(:,1));
  eigK_ = V*K(:,1);
  K_ = V*diag(sort(eigK,'ascend'))*conj(V')/npts;
  Ks = V*diag(sqrt(eigK))*conj(V')/npts;
end

b = 0.5; shape = 0.8;
% integrand = @(x) kernelSmooth(x,shape,b);
integrand = @(x) prod(2*b*((cos(2*pi*x)-b))./(1 + b^2 - 2*b*cos(2*pi*x)), 2);
%bernPoly = @(x)(-x.*(1-x) + 1/6); integrand = @(x) prod((1 + 0.9*bernPoly(x)), 2);

inputArgs = {'dim',2, 'absTol',1E-2, 'order',1, ....
      'stopAtTol',true, 'stopCriterion','MLE', 'arbMean',true };
inputArgs{end+1} = 'f'; inputArgs{end+1} = integrand;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = 'SmoothKernel';

obj=cubBayesLattice_g(inputArgs{:});
nRep = 100;
muhatVec(nRep) = 0;
%outVec(nRep) = 0;
for i=1:nRep
  [muhatVec(i),outVec(i)]=compInteg(obj);
end
median(muhatVec)


%[K,kvec,k0] = kernelFunBern(xlat,0.1,0.2);
[K,kvec,k0] = CovarMatrix(xlat_,shape,b);
figure; plot(xlat_, K(:,1)); axis tight
figure; scatter(xlat, K(:,npts/4), 10); axis tight

% check eigen values
ev = eig(K);
figure; loglog(abs(ev)); axis tight
ev2 = real(fft(K(:,1)));
figure; loglog(abs(sort(ev)./sort(ev2)))
if abs(sum(ev - sort(ev2))) > 1E-8
  warning('Fast transform: eigen values do not match');
end

% check positive definiteness
for ind=1:100
  c = rand(npts,1);
  if c'*K*c <= 0
    warning('Not positive definite');
  end
end

% y = bernPoly(xlat_);
% b = K*y;
% y_ = ifft(real(fft(b))./real(fft(K(:,1))));
% check matrix inverse
for ind=1:100
  c = rand(npts,1);
  b = K*c;
  c_ = K\b;
  if abs(sum(c-c_)) > 1E-8
    warning('Inverse: too big error !');
  end
end

% check the symmetric circulant property
K_ = toeplitz(K(:,1), K(1,:)');

for m=1:npts
  dg = diag(K,m-1)';
  if abs(sum(dg(1)-dg)) > 1E-8
    fprintf('%dth diagonal is wrong', m)
  end
end

plotObjectiveFunc(dim, z)

fprintf('done')
end

function plot_kernel()
  [X, Y] = meshgrid(l);

  for sh = shape_param
    for b = bvec
      hFig = figure('visible','on');
      set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
      Z = kernelSmooth([X(:) Y(:)], sh, b);
      meshc(X, Y, reshape(Z, [n, n]));

      title(sprintf('$b=%0.2f, \\theta$=%0.2f', b, sh), 'Interpreter','latex')
      figSavePathName = sprintf('fourier_kernel b_%dby100 shape_%dby100.png', 100*b, 100*sh);
      %saveas(hFig, figSavePathName)
    end
  end
end

function C1 = kernelSmooth(x,shape,b)
  tempa = 2*pi*x;
  % tempb = b*cos(tempa);
  % C1 = 1 + shape*(2*abs(tempb-1))./(1 + b^2 - 2*tempb);
  % C1 = 1 + shape*((b*b-1))./(1 + b^2 - 2*tempb);
  % C1 = ((b*b-1))./(1 + b^2 - 2*tempb) + 1;
  % C1 = 1 + shape*2*(b*(b-cos(tempa)))./(1 + b^2 - 2*b*cos(tempa));
  C1 = 1 + shape*2*b*((cos(tempa)-b))./(1 + b^2 - 2*b*cos(tempa));

end

function [K,kvec,k0] = CovarMatrix(x,shape,b)
[nx,d] = size(x);
k0 = 1;
kvec = ones(nx,1);
K = ones(nx);
for k = 1:d
  % tempa = abs(bsxfun(@minus,x(:,k),x(:,k)') );  % Lattice
  tempa = mod(bsxfun(@minus,x(:,k),x(:,k)'), 1);  % Lattice
  tempc = kernelSmooth(tempa,shape,b);
  K = K.*tempc;
end

end

function minTheta = plotObjectiveFunc(dim ,z)

obj.mvec = 8:16;
obj.dim = dim;
obj.ff = @(x) keisterFunc(x,obj.dim,1/sqrt(2)); % a=1/sqrt(2)
obj.fName='Keister';
obj.arbMean = true;
obj.order = 2;
obj.ptransform = 'Baker';
obj.figSavePath = '';


numM = length(obj.mvec);
n = 2.^obj.mvec(end);

%% plot ObjectiveFunction
lnTheta = -7:7;
% build filename with path to store the plot
% plotFileName = sprintf('%s%s_Cost_d%d_r%d_%s_%s.png',...
%   obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform, obj.stopCriterion);

costMLE = zeros(numM,numel(lnTheta));
tstart = tic;

% loop over all the m values
for iter = 1:numM
  nii = 2^obj.mvec(iter);
  
  xlat_ = mod(bsxfun(@times,(0:1/nii:1-1/nii)',z),1); % unshifted in direct order
  fx = obj.ff(xlat_);  % Note: periodization transform already applied
  ftilde_iter = fft(fx);
  br_xun = xlat_;
  
  tic
  %par
  for k=1:numel(lnTheta)
    [costMLE(iter,k)] = ObjectiveFunction(exp(lnTheta(k)),...
      br_xun,ftilde_iter);
  end
  toc
end

toc(tstart)

hFigCost = figure();

% semilogx
semilogx(exp(lnTheta),real(costMLE));
set(hFigCost, 'units', 'inches', 'Position', [1 1 10 7])

xlabel('Shape param, \(\theta\)')
ylabel('MLE Cost, \( \log \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
axis tight;
mType = '\(m \neq 0\)'; % arb mean

title(sprintf('%s d=%d r=%d Tx=%s %s', obj.fName, obj.dim, obj.order, obj.ptransform, mType));
[minVal,Index] = min(real(costMLE),[],2);

% mark the min theta values found using fminbnd
minTheta = exp(lnTheta(Index));
hold on;
semilogx(minTheta, minVal, '.');
if exist('aMLEAll', 'var')
  semilogx(obj.aMLEAll, obj.lossMLEAll, '+');
end
temp = string(obj.mvec);
temp = strcat('\(2^{',temp,'}\)');
temp(end+1) = '\(\theta_{min_{true}}\)';
if exist('aMLEAll', 'var')
  temp(end+1) = '\(\theta_{min_{est}}\)';
end
legend(temp,'location','best'); axis tight
plotFileName='MLE_onjective_func.png';
saveas(hFigCost, plotFileName)
end

function [loss,Lambda] = ObjectiveFunction(a,xun,ftilde)

n = length(ftilde);

b = 0.1;
C1 = prod(kernelSmooth(xun,a,b),2);
Lambda = abs(real(fft(C1)));
ftilde = abs(ftilde);  % remove any negative values

% compute RKHSnorm
temp = (abs(ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0))) ;

% compute loss
% default: MLE
temp_1 = sum(temp(2:end));

% ignore all zero eigenvalues
loss = sum(log(Lambda(Lambda~=0))) + n*log(temp_1);

end

function [K] = MaternkernelFun(x,shape)
[nx,d] = size(x);
K = ones(nx);
for k = 1:d
  tempa = shape*mod(bsxfun(@minus,x(:,k),x(:,k)'), 1);  % Lattice
  %tempa = shape*abs(bsxfun(@minus,x(:,k),x(:,k)'));
  K = K.*exp(-tempa).*(1 + tempa);
end

end
