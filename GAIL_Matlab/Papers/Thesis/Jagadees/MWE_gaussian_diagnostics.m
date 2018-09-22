%
% Minimum working example to test Gaussian diagnostics idea
%
function MWE_gaussian_diagnostics()
format short
close all
fName = 'ExpCos';
ptransform = 'C1sin';

npts = 2^10;  % max 14
dim = 1;
shift = rand(1,dim);
z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
            151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
            85729, 14597, 94813, 422013, 484367]; %generator      
z = z(1:dim);

xlat_ = mod(bsxfun(@times,(0:1/npts:1-1/npts)',z),1); % unshifted in direct order
xlat_ = mod(bsxfun(@plus,xlat_,shift),1);  % shifted in direct order

integrand = @(x) exp(sum(cos(2*pi*x), 2));
%integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8

integrand_p = cubBayesLattice_g.doPeriodTx(integrand, ptransform);

y = integrand_p(xlat_);
ftilde = fft(y);
ftilde(1) = 0;  % ((f - m_\MLE I)
if dim==1
figure; scatter(xlat_, y, 10)
end

bVec = [0.3]; %[0.3:0.1:0.9];  %
for sh=[0.001 0.01 0.1 0.5 0.9 5]
  for b=bVec
    C1 = kernel(xlat_,sh,b);
    lambda = fft(C1);
    
    w_ftilde = real(ifft(ftilde./sqrt(real(lambda))));
    hFigNormplot = figure();
    set(hFigNormplot,'defaultaxesfontsize',16, ...
      'defaulttextfontsize',16, ... %make font larger
      'defaultLineLineWidth',0.75, 'defaultLineMarkerSize',8)
    normplot(w_ftilde)
    
    % Shapiro-Wilk test
    [H, pValue, W] = swtest(w_ftilde);
    Hval='false';
    if H==true
      Hval='true';
    end
    fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3f, W=%1.3f\n', Hval, pValue, W);
    
    title(sprintf('%s n=%d Tx=%s b=%1.2f theta=%1.2f', ...
      fName, npts, ptransform, b, sh))
  end
end

end

function C1 = kernel(xun,whKer,a,b)

if strcmp(whKer,'Bern')
  b_order = 2;
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

C1 = prod(1 + (a)*constMult*kernelFunc(xun),2);  %

end

