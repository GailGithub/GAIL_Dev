%Run exact sampling for comparison 
%% QE_European call option
% %InitializeWorkspaceDisplay %initialize the workspace and the display parameters
clearvars
T=5;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T; 
% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties
initPrice = 100;
interest = 0.0;
inp.assetParam.initPrice = initPrice; %initial stock price
inp.assetParam.interest = interest; %risk-free interest rate
inp.assetParam.volatility = 0.3;
inp.assetParam.Vinst = 0.04;%0.36; %0.04; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
%inp.assetParam.nu = 0.5;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE';
inp.payoffParam.putCallType = {'call'};
%inp.priceParam.cubMethod = 'Sobol';
inp.priceParam.cubMethod = 'IID_MC';

%%
% To generate some discounted option payoffs to add some more properties
Strike = 90;
inp.payoffParam.strike =Strike; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.005; %one penny on the dollar relative tolerance
% ourGBMCallPrice = optPrice(inp)
% [GBMCallPrice_withoutIS, out] = genOptPrice(ourGBMCallPrice)
% return
Ntime = T/0.1; 
NSim = 1e3;
nu = [0,10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),0.05,10^(-1),0.15,0.2];
inp.assetParam.nu =0.5;% 0.1;%0.2;
%ExactSampling_Heston(S0,V0,K,r,T,kappa,theta,nu,rho,Nt,NSim)
% ExactSampling_Heston(initPrice,inp.assetParam.Vinst, Strike,interest,T,...
%     inp.assetParam.kappa,inp.assetParam.Vlong,inp.assetParam.nu,...
%     inp.assetParam.rho,Ntime,NSim)
% return
[ExactSamplingPrice_BroadieKaya,time] = HestonFullSampling(initPrice, Strike,interest,T,...
    inp.assetParam.kappa,inp.assetParam.Vlong,inp.assetParam.nu,...
    inp.assetParam.rho,inp.assetParam.Vinst,NSim,Ntime)
% Nt = T/delta_t;
% ExactSamplingPrice_BroadieKaya = ExactSampling_Heston(initPrice,inp.assetParam.Vinst,...
%     Strike,interest,T,inp.assetParam.kappa,inp.assetParam.Vlong,inp.assetParam.nu,...
%     inp.assetParam.rho,Nt,1e6)





