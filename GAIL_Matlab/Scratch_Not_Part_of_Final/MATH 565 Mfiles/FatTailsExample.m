% Fat tails
r=0.03; %interest rate
S0=100; %initial stock price
sig=0.5; %volatility
T=1; %expiry time
K=100; %strike price
n=1e6; %sample size

%% Normal random variables
disp('Normal (Gaussian)')
x=randn(n,1); %normal random
ST=exp(sig*sqrt(T)*x); %with std. dev. sig
ST=ST*(S0*exp(r*T)/mean(ST)); 
  %normalized to get expected value of ST right
callpay=max(ST-K,0)*exp(-r*T);
putpay=max(K-ST,0)*exp(-r*T);
callprice=mean(callpay);
putprice=mean(putpay);
callerr=1.96*std(callpay)/sqrt(n);
puterr=1.96*std(putpay)/sqrt(n);
disp(['call price = $' num2str(callprice,'%0.4f') char(177) num2str(callerr,'%0.4f')])
disp(['put price = $' num2str(putprice,'%0.4f') char(177) num2str(puterr,'%0.4f')])
disp(' ');

%% Uniform random variables
disp('Uniform')
x=sqrt(12)*(rand(n,1)-0.5); %uniform random
ST=exp(sig*sqrt(T)*x); %with std. dev. sig
ST=ST*(S0*exp(r*T)/mean(ST)); 
  %normalized to get expected value of ST right
callpay=max(ST-K,0)*exp(-r*T);
putpay=max(K-ST,0)*exp(-r*T);
callprice=mean(callpay);
putprice=mean(putpay);
callerr=1.96*std(callpay)/sqrt(n);
puterr=1.96*std(putpay)/sqrt(n);
disp(['call price = $' num2str(callprice,'%0.4f') char(177) num2str(callerr,'%0.4f')])
disp(['put price = $' num2str(putprice,'%0.4f') char(177) num2str(puterr,'%0.4f')])
disp(' ');

%% Binomial random variables
disp('Binomial')
x=2*(ceil(2*rand(n,1))-1); %Binomial random
ST=exp(sig*sqrt(T)*x); %with std. dev. sig
ST=ST*(S0*exp(r*T)/mean(ST)); 
  %normalized to get expected value of ST right
callpay=max(ST-K,0)*exp(-r*T);
putpay=max(K-ST,0)*exp(-r*T);
callprice=mean(callpay);
putprice=mean(putpay);
callerr=1.96*std(callpay)/sqrt(n);
puterr=1.96*std(putpay)/sqrt(n);
disp(['call price = $' num2str(callprice,'%0.4f') char(177) num2str(callerr,'%0.4f')])
disp(['put price = $' num2str(putprice,'%0.4f') char(177) num2str(puterr,'%0.4f')])
disp(' ');

%% T distribution random variables
deg=6;
disp('t-Distribution')
disp(['   with ' num2str(deg) ' degrees of freedom'])
x=sqrt((deg-2)/deg)*trnd(deg,n,1); %t distribution random
ST=exp(sig*sqrt(T)*x); %with std. dev. sig
ST=ST*(S0*exp(r*T)/mean(ST)); 
  %normalized to get expected value of ST right
callpay=max(ST-K,0)*exp(-r*T);
putpay=max(K-ST,0)*exp(-r*T);
callprice=mean(callpay);
putprice=mean(putpay);
callerr=1.96*std(callpay)/sqrt(n);
puterr=1.96*std(putpay)/sqrt(n);
disp(['call price = $' num2str(callprice,'%0.4f') char(177) num2str(callerr,'%0.4f')])
disp(['put price = $' num2str(putprice,'%0.4f') char(177) num2str(puterr,'%0.4f')])
disp(' ');

%% Variance gamma random variables
beta=1;
disp('Variance Gamma')
disp(['   with beta = ' num2str(beta)])
x=sqrt(gamrnd(T/beta,beta,n,1)).*randn(n,1); %variance gamma random
ST=exp(sig*x); %with std. dev. sig
ST=ST*(S0*exp(r*T)/mean(ST)); 
  %normalized to get expected value of ST right
callpay=max(ST-K,0)*exp(-r*T);
putpay=max(K-ST,0)*exp(-r*T);
callprice=mean(callpay);
putprice=mean(putpay);
callerr=1.96*std(callpay)/sqrt(n);
puterr=1.96*std(putpay)/sqrt(n);
disp(['call price = $' num2str(callprice,'%0.4f') char(177) num2str(callerr,'%0.4f')])
disp(['put price = $' num2str(putprice,'%0.4f') char(177) num2str(puterr,'%0.4f')])
disp(' ');

