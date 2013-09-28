close all, clear all, format short
d=5;
p=sobolset(d); %vanilla Sobol
psc=scramble(p,'MatousekAffineOwen'); %scrambled Sobol
n=256; %sample size
hx=1/32; hy=1/8; %sizes of the boxes
subplot(1,2,1); plot(p(1:n,1),p(1:n,2),'.','markersize',10)
set(gca,'Xtick',0:hx:1,'Ytick',0:hy:1); grid on
subplot(1,2,2); plot(psc(1:n,1),psc(1:n,2),'.','markersize',10)
set(gca,'Xtick',0:hx:1,'Ytick',0:hy:1); grid on
disp([p(1:5,1:2) psc(1:5,1:2)])

%% LHS
n=1e5;
d=5;
tic
xlhs=zeros(n,d);
for j=1:d
    xlhs(:,j)=2*(randperm(n)'-rand(n,1))/n-1; %create Latin Hypercuble sample
end
toc

%Or try it this way
tic; xMATlhs=lhsdesign(n,d,'criterion','none'); toc