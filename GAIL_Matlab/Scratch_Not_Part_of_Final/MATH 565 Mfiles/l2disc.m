function disc=l2disc(x)
[n,d]=size(x);

kmat=ones(n,n);
intkvec=ones(n,1);
c=1;
for j=1:d
  xmat=repmat(x(:,j),1,n);
  kmat=kmat.*(2 - max(xmat,xmat'));
  intkvec=intkvec.*((3-x(:,j).*x(:,j))/2);
  c=c*4/3;
end
disc=c-2*mean(intkvec)+mean(mean(kmat));
disc=sqrt(disc);