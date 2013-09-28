%See good lattice generators
clear all
n=50;
d=7;

h2=[1:n-1]; 
ngen=length(h2);
yj=repmat([0:n-1]',1,ngen);
kvec=ones(n,ngen);
for j=1:d;
    xj=yj/n;
    kvec=kvec.*(1+(1/6 + xj.*(-1 + xj)));
    yj=mod(yj.*repmat(h2,n,1),n);
end
disc2=-1+mean(kvec,1);
[best,wh]=min(disc2)
h2best=h2(wh)
