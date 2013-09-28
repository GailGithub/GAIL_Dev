%Simple pseudo-random number generator
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18)
format compact; clear all; close all
M=11;
a=6; nu=[4/11 2/11; 4/11-1/5 2/11+2/5]; plane=[0 1 0 0.5; 0 1 0.5 1]; npl=2;
%a=7; nu=[2/11 3/11; 2/11+1/10 3/11+3/10]; plane=[0 1 1/3 0; 0 1 2/3 1/3; 0 1 1 2/3]; npl=3;
m0=3;
for i=1:M-1
    m(i)=mod(a*m0,M);
    m0=m(i);
end
x=m/M;
meanx=mean(x)
varx=var(x)
xlag=x([2:M-1 1]);
ans=cov(x,xlag);
covxxlag=ans(1,2)
h=plot(x,xlag,'s'); set(h,'linewidth',6)
axis equal; axis([0 1 0 1])
hold on
xlabel('{\it x_{i}}');
ylabel('{\it x_{i+1}}');
eval(['print -depsc SmallLCGExM' int2str(M) 'a' int2str(a) '.eps'])

figure;
h=plot(x,xlag,'s'); set(h,'linewidth',6)
axis equal; axis([0 1 0 1])
hold on
xlabel('{\it x_{i}}');
ylabel('{\it x_{i+1}}');
eval(['print -depsc SmallLCGExM' int2str(M) 'a' int2str(a) '.eps'])
h=plot(nu(:,1),nu(:,2),'r-');  set(h,'linewidth',2)
for i=1:npl %draw planes through points
  h=plot(plane(i,[1 2]),plane(i,[3 4]),'k--');  set(h,'linewidth',2);
end
eval(['print -depsc SmallLCGExM' int2str(M) 'a' int2str(a) 'planes.eps'])