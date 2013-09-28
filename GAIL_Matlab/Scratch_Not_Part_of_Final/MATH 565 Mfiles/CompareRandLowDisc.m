%% Compare low discrepancy sequences
format compact, format long
d=2; %dimension
n=2^5; %sample size
%f=@(x) exp(sum(x,2)); %define function to be integrated
%muexact=(exp(1)-1)^d
f=@(x) 1+10*cos(sum(2*pi*x,2)); %define function to be integrated
muexact=1;
xrand=rand(n,d); %uniform random samples
xsob=net(scramble(sobolset(d),'MatousekAffineOwen'),n); %scrambled Sobol' sequence
xlat=mod(kuorank1(d,n)+repmat(rand(1,d),n,1),1); %shifted lattice
yrand=f(xrand);
ysob=f(xsob);
ylat=f(xlat);
yrandbar=mean(yrand) %answer for random points
ysobbar=mean(ysob) %answer for Sobol' points
ylatbar=mean(ylat) %answer for lattice points