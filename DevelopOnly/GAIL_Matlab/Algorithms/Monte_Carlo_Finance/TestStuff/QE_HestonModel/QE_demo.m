% This file is a demo for QE method compared with ExactSampling for Heston
% model.
% Latest updated date: 10/12/2016
% Author: Xiaoyang Zhao

%% QE_European call option
% %InitializeWorkspaceDisplay %initialize the workspace and the display parameters
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
inp.assetParam.Vinst = 0.09; %0.04; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
inp.assetParam.nu = 0; %0.8;
inp.assetParam.rho = 0; %-0.3;
inp.assetParam.pathType = 'GBM';
%inp.priceParam.cubMethod = 'Sobol';
inp.priceParam.cubMethod = 'IID_MC';

%%
% To generate some discounted option payoffs to add some more properties
Strike = 120;
inp.payoffParam.strike =Strike; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance
ourGBMCallPrice = optPrice(inp);
[GBMCallPrice, out] = genOptPrice(ourGBMCallPrice); %the option price

inp.assetParam.pathType = 'QE_m';
inp.assetParam.meanShift =0.49;
ourQECallPrice = optPrice(inp) %construct an optPrice object 
%genOptPayoffs(ourQECallPrice,1);
%return
[QECallPrice, out] = genOptPrice(ourQECallPrice) %the option price
inp.priceParam.cubMethod = 'Sobol'
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T; 
ourQECallPrice = optPrice(inp) %construct an optPrice object 
%genOptPayoffs(ourQECallPrice,1);
%return
[QECallPrice, out] = genOptPrice(ourQECallPrice) %the option price
return
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
% Ntime = numel(inp.timeDim.timeVector)-1;
tic
Ntime = T/0.25; 
NSim = 1e6;
[a,b] = MC_QE(initPrice,interest,0,T,inp.assetParam.Vinst,inp.assetParam.Vlong,...
    inp.assetParam.kappa,inp.assetParam.nu,inp.assetParam.rho,Ntime,NSim,1);
PT = a(:,Ntime + 1);
PT = max(PT-Strike,0);
PP = mean(PT);
QEprice_Kienitz = PP*exp(-inp.assetParam.interest*T)
toc
return
ExactSamplingPrice_BroadieKaya = HestonFullSampling(initPrice, Strike,interest,T,...
    inp.assetParam.kappa,inp.assetParam.Vlong,inp.assetParam.nu,...
    inp.assetParam.rho,inp.assetParam.Vinst,NSim,Ntime)

return
% Nt = T/delta_t;
% ExactSamplingPrice_BroadieKaya = ExactSampling_Heston(initPrice,inp.assetParam.Vinst,...
%     Strike,interest,T,inp.assetParam.kappa,inp.assetParam.Vlong,inp.assetParam.nu,...
%     inp.assetParam.rho,Nt,1e6)

inp.assetParam.pathType = 'QE_m';
ourQEmCallPrice = optPrice(inp) %construct an optPrice object 
%genOptPayoffs(ourQECallPrice,1);
%return
[QEmCallPrice, out] = genOptPrice(ourQEmCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
Ntime = T/delta_t; 
[a,b] = MC_QE_m(initPrice,interest,0,T,inp.assetParam.Vinst,inp.assetParam.Vlong,...
    inp.assetParam.kappa,inp.assetParam.nu,inp.assetParam.rho,Ntime,1e6,1);
PT = a(:,Ntime + 1);
PT = max(PT-Strike,0);
PP = mean(PT);
QEmprice_Kienitz = PP*exp(-inp.assetParam.interest*T)
return

%% QE_European call option_Strike = 70
%%
%InitializeWorkspaceDisplay %initialize the workspace and the display parameters
% inp.timeDim.startTime = 0;
% inp.timeDim.endTime = 5;
% inp.timeDim.nSteps  = 5;
inp.timeDim.timeVector = 0:1:5; 

% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties

inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0; %risk-free interest rate
inp.assetParam.volatility = 0.3;
inp.assetParam.Vinst = 0.09; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
inp.assetParam.nu = 0;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =70; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,nu,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE(100,0,0,5,0.09,0.09,0.0000001,0,-0.3,5,100000,1);
 Strike = 70;
 PT = a(:,6);
 PT = max(PT-Strike,0);
 PP = mean(PT)
%  PP = PP*exp(-0.01*3)

%% QE_European call option_Strike = 100 

%%
% %InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.endTime = 5; 
inp.timeDim.nSteps  = 15;

% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties

inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0; %risk-free interest rate
inp.assetParam.volatility = 0.3;
inp.assetParam.Vinst = 0.09; 
inp.assetParam.Vlong = 0.09;
inp.assetParam.kappa = 1;
inp.assetParam.nu = 0.0000001;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =100; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,nu,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE(100,0,0,5,0.09,0.09,1,0,-0.3,15,400000,1);
 Strike = 100;
 PT = a(:,16);
 PT = max(PT-Strike,0);
 PP = mean(PT)

%% QE_European call option_Strike = 105 

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
inp.assetParam.nu = 0;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =105; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,nu,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE(100,0,0,5,0.09,0.09,1,0,-0.3,15,700000,1);
 Strike = 105;
 PT = a(:,16);
 PT = max(PT-Strike,0);
 PP = mean(PT)

 %% QE_European call option_Strike = 140

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
inp.assetParam.nu = 0;
inp.assetParam.rho = -0.3;
inp.assetParam.pathType = 'QE';

%%
% To generate some discounted option payoffs to add some more properties

inp.payoffParam.strike =140; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance
ourCallPrice = optPrice(inp) %construct an optPrice object 

[CallPrice, out] = genOptPrice(ourCallPrice) %the option price
% Calculate option price by provided codes
%  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,nu,rho,NTime,NSim,NBatches)
 [a,b]=MC_QE(100,0,0,5,0.09,0.09,1,0,-0.3,15,1000000,1);
 Strike = 140;
 PT = a(:,16);
 PT = max(PT-Strike,0);
 PP = mean(PT)




