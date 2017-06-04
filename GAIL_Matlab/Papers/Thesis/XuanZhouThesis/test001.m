lambda=1;
kernel=@(x,t) exp(-lambda*(bsxfun(@minus,x,t')).^2);
n=20;
xnode=(0:(n-1))'/(n-1);
K=kernel(xnode,xnode);
[V,sig,~]=svd(K,0);
dsig=diag(sig);
m=find(dsig>1e-10,1,'last');
%B=V(:,1:m).*repmat(1./sqrt(dsig(1:m)'),n,1);
B = bsxfun(@times,V(:,1:m),1./sqrt(dsig(1:m)'));
nrep=1000;
ntest=100;
fval=zeros(ntest,nrep);
xtest=rand(ntest,1);
for i=1:nrep
   coef=B*randn(m,1);
   frand=@(x) kernel(x,xnode)*coef;
   fval(:,i)=frand(xtest);
end
obscov=cov(fval');
Ktest=kernel(xtest,xtest);
err=norm(Ktest-obscov)/norm(Ktest);
display(err);
