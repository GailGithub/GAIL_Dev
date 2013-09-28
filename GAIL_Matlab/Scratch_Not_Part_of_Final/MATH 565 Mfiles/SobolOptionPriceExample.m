% This example prices an option whose payoff is
%  max(max(S_1(T),S_2(T))-K,0)
clear all, close all
format compact
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font labels large enough

S10=100; %initial price of stock 1
S20=100; %initial price of stock 2
sig1=0.5; %volatility of stock 1
sig2=0.7; %volatility of stock 2
rho=0.6; %correlation of the two Brownian motions
r=0; %interest rate
T=1; %time to expiry
K=100; %strike price
d=2;

%% Simple IID Monte Carlo
ntot=2^20; %total number of samples
tic
xnorm=randn(ntot,d); %normal random numbers
brown=[xnorm(:,1) ... %correlated Brownian motion
    rho*xnorm(:,1)+sqrt(1-rho^2)*xnorm(:,2)];
%two stock prices
S1T=S10*exp((r-sig1^2/2)*T + sig1*sqrt(T)*brown(:,1));
S2T=S10*exp((r-sig2^2/2)*T + sig2*sqrt(T)*brown(:,2));
discpayoff=max(max(S1T,S2T)-K,0)*exp(-r*T); %payoff
priceiid=mean(discpayoff); %option price
erriid=2.58*std(discpayoff)/sqrt(ntot); %error estimate
timeiid=toc;
disp(['IID price = $' num2str(priceiid) ' +/- $' ...
    num2str(erriid) ' in ' num2str(timeiid) ' seconds'])
nplot=min(2^9,ntot);
figure;
plot(xnorm(1:nplot,1),xnorm(1:nplot,2),'b.')
axis([-3 3 -3 3]); axis('square')
set(gca,'Xtick',-3:3,'Ytick',-3:3)

%% Scrambled Sobol' with IID replications
%Central limit theorem applies
ntot=2^13; %total number of samples
nrep=2^4; %number of replications of random Sobol'
nSob=ntot/nrep; %number of Sobol' points per replication
tic
priceoneSob=zeros(nrep,1);
for i=1:nrep
    scsobol=qrandstream(scramble(sobolset(d),'MatousekAffineOwen'));
    xnorm=norminv(rand(scsobol,nSob,d)); %inverse normal transformation
    brown=[xnorm(:,1) ...  %correlated Brownian motion
        rho*xnorm(:,1)+sqrt(1-rho^2)*xnorm(:,2)];
    %two stock prices
    S1T=S10*exp((r-sig1^2/2)*T + sig1*sqrt(T)*brown(:,1));
    S2T=S10*exp((r-sig2^2/2)*T + sig2*sqrt(T)*brown(:,2));
    discpayoff=max(max(S1T,S2T)-K,0)*exp(-r*T); %payoff
    priceoneSob(i)=mean(discpayoff); %option price for one replication
end
priceSobIID=mean(priceoneSob);
errSobIID=2.58*std(priceoneSob)/sqrt(nrep); %error estimate
timeSobIID=toc;
disp(['Sobol IID rep price = $' num2str(priceSobIID) ' +/- $' ...
        num2str(errSobIID) ' in ' num2str(timeSobIID) ' seconds'])
figure;
nplot=min(2^9,nSob);
plot(xnorm(1:nplot,1),xnorm(1:nplot,2),'b.')
axis([-3 3 -3 3]); axis('square')
set(gca,'Xtick',-3:3,'Ytick',-3:3)
    
%Scrambled Sobol' with internal replications
%No theory, just heuristic
tic
scsobol=qrandstream(scramble(sobolset(d),'MatousekAffineOwen'));
xnorm=norminv(rand(scsobol,ntot,d)); %inverse normal transformation
brown=[xnorm(:,1) ...  %correlated Brownian motion
    rho*xnorm(:,1)+sqrt(1-rho^2)*xnorm(:,2)];
%two stock prices
S1T=S10*exp((r-sig1^2/2)*T + sig1*sqrt(T)*brown(:,1));
S2T=S10*exp((r-sig2^2/2)*T + sig2*sqrt(T)*brown(:,2));
discpayoff=max(max(S1T,S2T)-K,0)*exp(-r*T); %payoff
priceoneSob=mean(reshape(discpayoff,nSob,nrep),1); %option price for one replication
priceSobint=mean(priceoneSob);
errSobint=2.58*std(priceoneSob)/sqrt(nrep); %error estimate
timeSobint=toc;
disp(['Sobol internal rep price = $' num2str(priceSobint) ' +/- $' ...
        num2str(errSobint) ' in ' num2str(timeSobint) ' seconds'])
figure;
plot(xnorm(1:nplot,1),xnorm(1:nplot,2),'b.')
axis([-3 3 -3 3]); axis('square')
set(gca,'Xtick',-3:3,'Ytick',-3:3)


