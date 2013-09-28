%Pictures of stock prices
%% Initialization
close all, clear all
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18)
tic; %start timer
N=2000; %number of stock paths
Nmax=50; %maxinum number of stock paths to plot
s0=100; %initial stock price
strike=80; %option strike price
r=0.05; %interest rate
sig=0.3; %volatility
d=12; %number of time steps
delt=1/d; %time step

%% Generate stock paths
x=randn(N,d); %generate normal random numbers
%generate a matrix containing N stock paths
smat=cumprod([s0*ones(N,1) exp((r-sig.*sig/2)*delt + x*sig*sqrt(delt))],2);
tvec=(0:d)*delt; %vector of times
plot(tvec,smat(1:min(N,Nmax),:),'-','linewidth',2) %plot the stock paths
xlabel('Time')
ylabel('Stock Price')
eval(['print -depsc MCStockN' int2str(N) '.eps'])

%% Compute payoff and option price
%European call option
orderedpay=(sort(smat(:,d+1))-strike)*exp(-r); %ordered discounted Asian call payoff if exercised
%Asian call option
%avgS=sum(smat(:,2:d+1),2)/(d+1); %average of stock prices
%orderedpay=(sort(avgS)-strike)*exp(-r); %ordered discounted Asian call payoff if exercised
payoff=max(orderedpay,0); %payoff of option
price=mean(payoff); %option price
err=1.96*std(payoff)/sqrt(N); %estimated error bound
disp(['Fair price of option = $' num2str(price,'%4.2f') ' +/- $' num2str(err,'%4.2f')])
prob=((1:N)-0.5)/N; %probabilities
priceprob=interp1(orderedpay,prob,price,'linear'); %interpolate to find where price is
figure; 
plot(payoff,((1:N)-0.5)/N,'-b',...
    [price,price],[0 priceprob],'k--',price,priceprob,'r.',...
    'linewidth',2,'markersize',30)
xlabel('Payoff');
ylabel('Probability');
text(price*1.1,0.1,'Fair Price')
axis([0 80 0 1])
eval(['print -depsc MCMeanCallN' int2str(N) '.eps'])
display(['Time for simulation and display = ' num2str(toc) ' seconds'])

