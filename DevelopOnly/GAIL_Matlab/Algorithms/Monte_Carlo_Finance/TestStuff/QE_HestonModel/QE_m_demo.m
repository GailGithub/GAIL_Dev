% Latest updated date: 4/7/2017
% Auther: Xiaoyang Zhao
%% QE with martingale correction_European call option_Strike = 70 
clearvars
T=5;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T; 
% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties
initPrice = 60;
interest = 0.0;
inp.assetParam.initPrice = initPrice; %initial stock price
inp.assetParam.interest = interest; %risk-free interest rate
inp.assetParam.volatility = 0.4;
inp.assetParam.Vinst = 0.5; 
inp.assetParam.Vlong = 0.16;
inp.assetParam.kappa = 1;
inp.assetParam.nu = 0.4;%1e-16;
inp.assetParam.rho = -0.3;
inp.priceParam.cubMethod = 'IID_MC';
% inp.priceParam.cubMethod = 'Sobol';
inp.payoffParam.putCallType = {'put'};
inp.payoffParam.optType = {'american'};
% inp.assetParam.pathType = 'GBM';

%%
% To generate some discounted option payoffs to add some more properties
Strike =20;
inp.payoffParam.strike =Strike; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance
% ourGBMCallPrice = optPrice(inp);
% [GBMCallPrice, out] = genOptPrice(ourGBMCallPrice); %the option price
inp.assetParam.meanShift = -0.8;%1;

inp.assetParam.pathType = 'QE';
ourQEmCallPrice = optPrice(inp) %construct an optPrice object 
%genOptPayoffs(ourQECallPrice,1);
%return
[QEmCallPrice, out] = genOptPrice(ourQEmCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
return
Ntime = T/delta_t; 
tic
[a,b] = MC_QE_m(initPrice,interest,0,T,inp.assetParam.Vinst,inp.assetParam.Vlong,...
    inp.assetParam.kappa,inp.assetParam.nu,inp.assetParam.rho,Ntime,1e6,1);
PT = a(:,Ntime + 1);
if strcmp(inp.payoffParam.putCallType,'call')
    PT = max(PT-Strike,0);
else
    PT = max(Strike-PT,0);
end

PP = mean(PT);
QEmprice_Kienitz = PP*exp(-inp.assetParam.interest*T)
elapsed = toc;
elapsed
return

%% QE with martingale correction_European call option_Strike = 100 

%%
% %InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.endTime = 5; 
inp.timeDim.nSteps  = 15;

% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties

inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0; %risk-free interest rate
inp.assetParam.volatility = 0.09;
inp.assetParam.Vinst = 0.09; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
inp.assetParam.epsilon = 0;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE_m';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =100; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.005; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE_m(100,0,0,5,0.09,0.09,1,0,-0.3,15,1e5,1);
 Strike = 100;
 PT = a(:,16);
 PT = max(PT-Strike,0);
 PP = mean(PT)

%% QE with martingale correction_European call option_Strike = 105 

%%
% %InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.endTime = 5; %time increments of 0.004 up to time 1
inp.timeDim.nSteps  = 15;

% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties

inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0; %risk-free interest rate
inp.assetParam.volatility = 0.09;
inp.assetParam.Vinst = 0.09; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
inp.assetParam.epsilon = 0;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE_m';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =105; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.005; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE_m(100,0,0,5,0.09,0.09,1,0,-0.3,15,100000,1);
 Strike = 105;
 PT = a(:,16);
 PT = max(PT-Strike,0);
 PP = mean(PT)

 %% QE with martingale correction_European call option_Strike = 140

%%
% %InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.endTime = 5; %time increments of 0.004 up to time 1
inp.timeDim.nSteps  = 15;

% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties

inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0; %risk-free interest rate
inp.assetParam.volatility = 0.09;
inp.assetParam.Vinst = 0.09; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
inp.assetParam.epsilon = 0;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE_m';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =140; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.005; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE_m(100,0,0,5,0.09,0.09,1,0,-0.3,15,100000,1);
 Strike = 140;
 PT = a(:,16);
 PT = max(PT-Strike,0);
 PP = mean(PT)
 
 