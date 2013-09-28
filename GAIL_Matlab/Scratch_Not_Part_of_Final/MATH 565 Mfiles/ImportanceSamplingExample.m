%Importance Sampling example
n=1e4;

g=@(x) x.*x;
gtilde=@(z) z/2;
mu=1/3

%Simple MC
tic
n=4e4;
x=rand(n,1);
y=g(x);
muhat=mean(y)
err=1.96*std(y)/sqrt(n)
toc

%Importance Sampling
tic
n=4e4;
z=sqrt(rand(n,1));
w=gtilde(z);
muhatimp=mean(w)
errimp=1.96*std(w)/sqrt(n)
toc