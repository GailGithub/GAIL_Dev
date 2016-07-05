function [K,kvec] = kernelFun(x,whKer)
x = gpuArray(x);
t = tic;
if nargin < 2
   whKer = 'sqExp';
end
[nx,d] = size(x);
K = ones(nx);
if strcmp(whKer,'sqExp')
   aa = 1;
   kvec = ones(nx,1)*(sqrt(pi)/(2*aa))^d;
   for k = 1:d;
      K = K.*exp(-(aa*bsxfun(@minus,x(:,k),x(:,k)')).^2);
      kvec = kvec.*(erf(aa*x(:,k)) + erf(aa*(1 - x(:,k))));
   end
elseif strcmp(whKer,'Mat1')
   aa = nx/100;
   %aa = 2;
   kvec = ones(nx,1)*((2/aa)^d);
   for k = 1:d;
      tempa = aa*abs(bsxfun(@minus,x(:,k),x(:,k)'));
      K = K.*exp(-tempa).*(1 + tempa);
      tempb = aa*x(:,k);
      tempc = aa - tempb;
      kvec = kvec.*(2 - exp(-tempc).*(1+tempc/2) ...
          - exp(-tempb).*(1+tempb/2));
   end
end
    toc(t)
end
