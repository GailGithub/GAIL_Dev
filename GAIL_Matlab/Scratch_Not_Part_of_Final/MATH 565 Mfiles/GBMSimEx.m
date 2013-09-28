N=1000000;
t=100;
x=rand(N,1); 
y=sqrt(2)*erfinv(2*x-1);
z=exp(sqrt(t)*y-t/2);

meanz=mean(z)