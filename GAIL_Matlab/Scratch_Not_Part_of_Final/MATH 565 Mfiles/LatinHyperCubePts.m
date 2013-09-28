%Latin Hypercube Sampling
close all, format compact, format short
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18)

n=50; %total number of points
m=5; %number of replications
k=n/m; %number of strata

%Simple Random
x=rand(n,2);
figure
h=plot(x(:,1),x(:,2),'rs'); set(h,'linewidth',6)
axis equal; axis([0 1 0 1]); 
xlabel('\it{x_1}')
ylabel('\it{x_2}')
set(gca,'xtick',[0:k]/k,'ytick',[0:k]/k,'xticklabel',[],'yticklabel',[])
grid on
eval(['print -depsc SimpRan' int2str(n) '.eps'])

%Latin Hypercube
xhyp=x;
for j=1:2;
    xtemp=reshape(x(:,j),k,m);
    if j>1; 
        for l=1:m; rpmat(:,l)=randperm(k)-1; end
    else
        rpmat=repmat([0:k-1]',1,m);
    end
    xtemp=(rpmat+xtemp)/k;
    xhyp(:,j)=reshape(xtemp,n,1);
end
figure
h=plot(xhyp(:,1),xhyp(:,2),'rs'); set(h,'linewidth',6)
axis equal; axis([0 1 0 1]); 
xlabel('\it{x_1}')
ylabel('\it{x_2}')
set(gca,'xtick',[0:k]/k,'ytick',[0:k]/k,'xticklabel',[],'yticklabel',[])
grid on
eval(['print -depsc SimpRan' int2str(n) '.eps'])

%Toy Example
fun=@(x,a)(1 + a.*(x(:,1)-1/2)).*(1+a.*(x(:,2)-1/2));
intf=1;
disp('      a       MC est   LHS est  MC err / LHS err')
for i=[-3:3];
    a=4.^i;
    MCint=mean(fun(x,a));
    LHint=mean(fun(xhyp,a));
    errratio= abs((1-MCint)/(1-LHint));
    disp([a MCint LHint errratio])
end

