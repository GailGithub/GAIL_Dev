function [vSol,out]=cubBayesLatticeCLASS(varargin)
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
tstart = tic; %start the clock
inp = gail.cubBayesLatticeParam(varargin{:});
out = gail.cubBayesLatticeOut(inp);

d = out.d;
z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
   151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
   85729, 14597, 94813, 422013, 484367]; %generator
z = z(1:d);
mmin = ceil(log2(out.nInit)); %minimum power of 2
mmax = ceil(log2(out.nMax)); %maximum power of 2
mvec = mmin:mmax;
numM = length(mvec);
if out.isShift
   shift = rand(1,d);
end
f = out.fff;
for ii = 1:numM
   m = mvec(ii);
   n = 2^m;
   
   %Update function values
   if ii == 1
      xun = mod(bsxfun(@times,(0:1/n:1-1/n)',z),1);
      if out.isShift
         x = mod(bsxfun(@plus,xun,shift),1);
      else
         x = xun;
      end
      fx = f(x);
      ftilde = fft(fx);
   else
      xunnew = mod(bsxfun(@times,(1/n:2/n:1-1/n)',z),1);
      if out.isShift
         xnew = mod(bsxfun(@plus,xunnew,shift),1);
      else
         xnew = xunnew;
      end
      temp = zeros(n,d);
      temp(1:2:n-1,:) = xun;
      temp(2:2:n,:) = xunnew;
      xun = temp;
      temp(1:2:n-1,:) = x;
      temp(2:2:n,:) = xnew;
      x = temp;
      fnew = f(xnew);
      fx = reshape([fx fnew]',n,1);
      ftilde = fft(fx);
   end
   
   %Compute MLE parameter
   lnaMLE = fminbnd(@(lntheta) ...
      MLEKernel(exp(lntheta),xun,ftilde,out), -5,5,optimset('TolX',1e-2));
   thetaMLE = exp(lnaMLE);
   [~,Ktilde,RKHSnorm] = MLEKernel(thetaMLE,xun,ftilde,out);
   
   %Check error criterion
   out.errBd = 2.58*sqrt(((1/n) - 1./Ktilde(1))*RKHSnorm);
   if ~numel(out.GPMean)
      muhat = ftilde(1)/n;
   else
      muhat = ftilde(1)/Ktilde(1);
   end
   vpm = out.solBdFun(muhat, out.errBd);
   tolPc = max(out.absTol,out.relTol*abs(vpm));
   den = tolPc(1) + tolPc(2);
   tolVal = abs(diff(vpm))/den;
   
   if tolVal <= 1
      out.tolVal = tolVal;
      vSol = (vpm(1) * tolPc(2) + vpm(2) * tolPc(1))/den;
      break
   end 
   
end
out.nSample = n;
out.mu = vSol;
out.time = toc(tstart);
end

function [loss,Ktilde,RKHSnorm] = MLEKernel(theta,xun,ftilde,out)
   K = out.kernel(xun,theta); %evaluate the kernel at the unshifted points
   Ktilde = real(fft(K)); %FFT of the values of the kernel
   temp = abs(ftilde).^2./Ktilde;
   if ~numel(out.GPMean)
      RKHSnorm = sum(temp(2:end))/numel(Ktilde);
   else
      RKHSnorm = mean(temp);
   end
   loss = mean(log(Ktilde)) + log(RKHSnorm);
end
