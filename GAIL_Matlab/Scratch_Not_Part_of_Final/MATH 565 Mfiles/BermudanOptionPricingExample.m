% Bermudan Option Pricing Example

%% Inititalize
clear all, close all
format compact

reset(RandStream.getGlobalStream)

n=1000;
nshow=10;
%rand(200,1);
S0=70;
K=100;
d=4;
T=1;
r=0.02;
sigma=0.7;
delt=T/d;
tvec=(0:d)*delt;

%% Generate stock price paths
Spath=zeros(n,d+1);
Spath(:,1)=S0;
Spath(:,2:d+1)=exp((r-sigma^2/2)*delt+sigma*sqrt(delt)*randn(n,d));
Spath=cumprod(Spath,2);
Spath(1:nshow,:)
payoff=max(K-Spath,0).*repmat(exp(-r*tvec),n,1);
payoff(1:nshow,:)

%% Compute exercise boundary
exerBdy=zeros(1,d+1);
exerBdy(d+1)=K
exerciseTime=repmat(d,n,1);
value=zeros(n,d+1); 
value(:,d+1)=payoff(:,d+1);
hold=zeros(n,d);
for j=d:-1:1
    whichInMoney=payoff(:,j)>0; %which paths in the money
    valueReg=value(whichInMoney,j+1);
    regMat=[ones(sum(whichInMoney),1) Spath(whichInMoney,j)];
    regCoef=regMat\valueReg
    hold(:,j)=[ones(n,1) Spath(:,j)]*regCoef;
    value(:,j)=max(payoff(:,j),hold(:,j));
    exerBdy(j)=(K*exp(-r*tvec(j))-regCoef(1))/(regCoef(2)+exp(-r*tvec(j)))
    whichExercise=(payoff(:,j)>hold(:,j))&(payoff(:,j)>0);
    exerciseTime(whichExercise)=j-1;
    [Spath(1:nshow,j) payoff(1:nshow,j) hold(1:nshow,j) value(1:nshow,j) exerciseTime(1:nshow)]
    keyboard
end
    
    
    
