function [muhat,out]=cubBayesLattice(f,d,absTol,relTol,order,arbMean)
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
   arbMean = true;
   if nargin < 5
      order = 2; %type of sampling, scrambled Sobol
      if nargin < 4
         relTol = 0;
         if nargin < 3
            absTol = 0.01; %dimension
            if nargin < 2
               d = 1; %number of samples
               if nargin < 1
                  f = @(x) x.^2; %function
               end
            end
         end
      end
   end
end

tstart = tic; %start the clock
z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
   151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
   85729, 14597, 94813, 422013, 484367]; %generator
z = z(1:d);
mmin = 10;
mmax = 23;
mvec = mmin:mmax;
numM = length(mvec);
shift = rand(1,d);
% ff = f; %no transformation
ff = @(x) f(1 - abs(1-2*x)); %folding transformation
for ii = 1:numM
   m = mvec(ii);
   n = 2^m;
   
   %Update function values
   if ii == 1
      xun = mod(bsxfun(@times,(0:1/n:1-1/n)',z),1);
      x = mod(bsxfun(@plus,xun,shift),1);
      fx = ff(x);
      ftilde = fft(fx);
   else
      xunnew = mod(bsxfun(@times,(1/n:2/n:1-1/n)',z),1);
      xnew = mod(bsxfun(@plus,xunnew,shift),1);
      temp = zeros(n,d);
      temp(1:2:n-1,:) = xun;
      temp(2:2:n,:) = xunnew;
      xun = temp;
      temp(1:2:n-1,:) = x;
      temp(2:2:n,:) = xnew;
      x = temp;
      fnew = ff(xnew);
      fx = reshape([fx fnew]',n,1);
      ftilde = fft(fx);
   end
   
   %Compute MLE parameter
   lnaMLE = fminbnd(@(lna) ...
      MLEKernel(exp(lna),xun,ftilde,order,arbMean), ...
      -5,5,optimset('TolX',1e-2));
   aMLE = exp(lnaMLE);
   [~,Ktilde,RKHSnorm] = MLEKernel(aMLE,xun,ftilde,order,arbMean);
   
   %Check error criterion
   out.ErrBd = 2.58*sqrt(((1/n) - 1./Ktilde(1))*RKHSnorm);
   if arbMean
      muhat = ftilde(1)/n;
   else
      muhat = ftilde(1)/Ktilde(1);
   end
   muminus = muhat - out.ErrBd;
   muplus = muhat + out.ErrBd;
   
   if 2*out.ErrBd <= ...
         max(absTol,relTol*abs(muminus)) + max(absTol,relTol*abs(muplus))      
      break
   end 
   
end
out.n = n;
out.time = toc(tstart);
end

function [loss,Ktilde,RKHSnorm] = MLEKernel(a,xun,ftilde,order,arbMean)
   if order == 2
      K = prod(1 + a*(-xun.*(1-xun) + 1/6),2);
   elseif order == 4
      K = prod(1 - a*((xun.*(1-xun)).^2 - 1/30),2);      
   end
   Ktilde = real(fft(K));
   temp = abs(ftilde).^2./Ktilde;
   if arbMean
      RKHSnorm = mean(temp(2:end));
   else
      RKHSnorm = mean(temp);
   end
   loss = mean(log(Ktilde)) + log(RKHSnorm);
%    a;
%    ftilde(1)/Ktilde(1);
end
