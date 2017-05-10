% Latest updated date: 4/7/2017
% Auther: Xiaoyang Zhao
%% QE with martingale correction_European call option_Strike = 70 
clearvars
T=5;
delta_t=0.5;
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
% inp.priceParam.cubMethod = 'IID_MC';
inp.priceParam.cubMethod = 'Sobol';
inp.payoffParam.putCallType = {'put'};
% inp.payoffParam.putCallType = {'call'};
inp.payoffParam.optType = {'american'};
% inp.assetParam.pathType = 'GBM';

%%
% To generate some discounted option payoffs to add some more properties
Strike = 20;
inp.payoffParam.strike =Strike; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.001; %one penny on the dollar relative tolerance
% ourGBMCallPrice = optPrice(inp);
% [GBMCallPrice, out] = genOptPrice(ourGBMCallPrice); %the option price

inp.assetParam.pathType = 'QE';
% inp.assetParam.meanShift = -0.6;%1;
ourQEmCallPrice = optPrice(inp) %construct an optPrice object 

avg = 1;
priceIS = 0;
timeIS = 0;
for i = 1:avg
    [QEmCallPrice, out] = genOptPrice(ourQEmCallPrice) %the option price
    priceIS = priceIS + QEmCallPrice;
    timeIS = timeIS + out.time;
end
priceIS = priceIS/avg;
timeIS = timeIS/avg;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanshifts = -0.1:-0.1:-1.2;
% meanshifts = 0.1:0.1:1.2;
avg = 10;
content = zeros(4,9);
Row = 1;
Col = 0;
for j = 1: length(meanshifts)
    inp.assetParam.meanShift = meanshifts(j);%1;
    ourQEmCallPrice = optPrice(inp) %construct an optPrice object 
    priceIS = 0;
    timeIS = 0;
    prices = zeros(avg,1);
    for i = 1:avg
        [QEmCallPrice, out] = genOptPrice(ourQEmCallPrice); %the option price
        priceIS = priceIS + QEmCallPrice;
        timeIS = timeIS + out.time;
        prices(i) = QEmCallPrice;
    end
    priceIS = priceIS/avg;
    timeIS = timeIS/avg;
    varianceIS = var(prices);
    if (mod(j-1,3)==0 && j~=1)
        Row = Row + 1;
        content(Row,1:3) = [priceIS,timeIS,varianceIS];
        Col = 1;
    else
        content(Row,1+3*Col:3*(Col+1)) = [priceIS,timeIS,varianceIS];
        Col = Col + 1;
    end
end
%% Generate LaTex code for a table of "call and put options",pathtype='QE',cubmethod='Sobol'
clear input;
% temp = NaN(1,6);
% call_QE = [IIDCallPrice_QE,temp];%,SobolCallPrice_QE,temp];
% put_QE = [IIDPutPrice_QE,temp];%,SobolPutPrice_QE,temp];
% content = [call_QE;Sobol_temp(1:2,:); put_QE;Sobol_temp(3:end,:)];
% content = [Sobol_temp(1:2,:);Sobol_temp(3:end,:)];
% content = [call_QE;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_QE;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'Price','Time','Variance','Price','Time','Variance','Price','Time','Variance'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'APw/IS','APw/IS','APw/IS','APw/IS'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat={'%.4f',1,'%.4f',1,'%3.2e',1,'%.4f',1,'%.4f',1,'%3.2e',1,'%.4f',1,'%.4f',1,'%3.2e',1};%,'%.4f',3,'%5.2fs',3,'%3.2E',3};
input.dataNanString = '-';
input.tableColumnAlignment = 'c';
% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'Sobol S0=60';
% LaTex table label:
input.tableLabel = 'MyTableLabel';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Switch to landscape table:
input.landscape = 0;
% call latexTable:
latex = latexTable(input);
%}

return

inp.assetParam.meanShift = 0;
ourQEmCallPrice = optPrice(inp)
price = 0;
time = 0;
prices = zeros(avg,1);
for i = 1:avg
    [QEmCallPrice, out] = genOptPrice(ourQEmCallPrice) %the option price
    price = price + QEmCallPrice;
    time = time + out.time;
    prices(i) = QEmCallPrice;
end
price = price/avg;
time = time/avg;
variance = var(prices);
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
 
 