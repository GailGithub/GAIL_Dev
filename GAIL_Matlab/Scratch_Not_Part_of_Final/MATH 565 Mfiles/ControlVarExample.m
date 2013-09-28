%% Control variates example

n=1000;

x=rand(n,1);
y=exp(x);
one=ones(n,1);
x2=x.*x;
meany=mean(y);
mean1=mean(one);
mean2=mean(x);
mean3=mean(x2);
X=[one-mean1 x-mean2 x2-mean3];
yvec=y-meany;
beta=X\yvec;

muhat=meany...
    -beta(1)*(mean1-1)...
    -beta(2)*(mean2-1/2)...
    -beta(3)*(mean3-1/3);