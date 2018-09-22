%
% Minimum working example to test if Maten kernel with Lattice points 
%    forms a Toeplitz Gram matrix
%
function MWE_Matern_Toeplitz()
format short
npts = 2^12;  % max 14
dim = 2;
shift = rand(1,dim);
z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
            151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
            85729, 14597, 94813, 422013, 484367]; %generator      
z = z(1:dim);
      
xlat_ = mod(bsxfun(@times,(0:1/npts:1-1/npts)',z),1); % unshifted in direct order
xlat = mod(bsxfun(@plus,xlat_,shift),1);  % shifted in direct order

%[K,kvec,k0] = kernelFunBern(xlat,0.1,0.2);
[K,kvec,k0] = kernelFunBern(xlat_,2,0.2);
figure; plot(xlat_, K(:,1)); axis tight
figure; plot(xlat, K(:,npts/4)); axis tight

K = kernelFun(xlat,0.1);

c = rand(npts,1);
b = K*c;
c_ = K\b;
sum(c-c_)
c'*K*c
K_ = toeplitz(K(:,1), K(1,:)');

for m=1:npts
  dg = diag(K,m-1)';
  if abs(sum(dg(1)-dg)) > 1E-8
    fprintf('%dth diagonal is wrong', m)
  end
end

fprintf('done')
end

function [K] = kernelFun(x,shape)
[nx,d] = size(x);
K = ones(nx);
for k = 1:d
  tempa = shape*mod(bsxfun(@minus,x(:,k),x(:,k)'), 1);  % Lattice
  %tempa = shape*abs(bsxfun(@minus,x(:,k),x(:,k)'));
  K = K.*exp(-tempa).*(1 + tempa);
end

end

function [K,kvec,k0] = kernelFunBern(x,shape,b)
[nx,d] = size(x);
k0 = 1;
kvec = ones(nx,1);
K = ones(nx);
for k = 1:d
  tempa = 2*pi*abs(bsxfun(@minus,x(:,k),x(:,k)') );  % Lattice
  tempb = b*cos(tempa);
  tempc = 1 + shape*(2*(tempb-1))./(1 + b^2 - 2*tempb);
  K = K.*tempc;
end

end
