%Test out nu fft
close all, clear all
format long
format compact


%% We test that our nufft gives us the same result as the Fourier transform
f=@(x) sin(x).*exp(cos(x));
n=2^15;
latticeseq_b2('init0');
xpts=latticeseq_b2(1,n);
ynu=f(xpts);
[Y,I]=sort(xpts)
yff=f(sort(xpts));
nfftr=nufft(2,ynu);
fftr=fft(yff)/n;
error=norm(nfftr-fftr);
plot(abs(nfftr-fftr))
figure
loglog(1:n,abs(nfftr),'r-',1:n,abs(fftr),'b-','linewidth',2)
hold on
loglog(1:n,abs(nfftr),'r-','linewidth',2)
loglog(1:n,abs(fftr),'b-','linewidth',2)
