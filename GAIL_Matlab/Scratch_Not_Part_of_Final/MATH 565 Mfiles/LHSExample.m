%Small Latin Hypercube Example
clear all
format short e
d=10;
f=@(x) sum(x,2)+...
    0.1*sum(x(:,1:d-1).*x(:,2:d),2);
exactint=d/2+0.1*(d-1)/4
n=1e3;
xrand=rand(n,d);
muMC=mean(f(xrand))
errMC=abs(exactint-muMC)

xlhs=lhsdesign(n,d);
mulhs=mean(f(xlhs))
errlhs=abs(exactint-mulhs)
