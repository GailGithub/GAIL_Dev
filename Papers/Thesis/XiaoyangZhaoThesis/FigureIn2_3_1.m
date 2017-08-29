% QE_European call option

%%InitializeWorkspaceDisplay %initialize the workspace and the display parameters
clearvars
T=2;%2;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T; 
% To generate an asset path modeled by a geometric Brownian motion we need
% to add some more properties
initPrice = 80;
interest = 0.0;
inp.assetParam.initPrice = initPrice; %initial stock price
inp.assetParam.interest = interest; %risk-free interest rate
inp.assetParam.volatility = 0.3;
inp.assetParam.Vinst = 0.36; %0.04; 
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
Strike = 100;
inp.payoffParam.strike =Strike; 

%% 
inp.priceParam.absTol = 0; %absolute tolerance
inp.priceParam.relTol = 0.01; %one penny on the dollar relative tolerance

Ntime = T/0.2; 
NSim = 1e6;
nu = [0,10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),0.05,10^(-1),0.2];
%nu = [0,10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),0.05,10^(-1),0.5];
ourQEprice = zeros(size(nu));
QEprice = zeros(size(nu));
for i = 1:length(nu)
    inp.assetParam.nu =nu(i);
    ourGBMCallPrice = optPrice(inp);
    [ourQEprice(i), out] = genOptPrice(ourGBMCallPrice);
    [pathS,pathV] = MC_QE(initPrice,interest,0,T,inp.assetParam.Vinst,inp.assetParam.Vlong,...
    inp.assetParam.kappa,inp.assetParam.nu,inp.assetParam.rho,Ntime,NSim,1);
    PT = pathS(:,Ntime + 1);
    PT = max(PT-Strike,0);
    PP = mean(PT);
    QEprice(i) = PP*exp(-inp.assetParam.interest*T);
end
dif = abs(ourQEprice-QEprice);
v0 = ourQEprice(1)*ones(size(nu));
QEzero = zeros(size(nu));
%}

%nu = log(nu);
%nu = exp(nu);
%QEprice = log(QEprice);
%%{
n=9;
figure
semilogx(nu(1:n),QEprice(1:n),'-o',nu(1:n),ourQEprice(1:n),'-*','LineWidth',1.5)
axis([0 1 -1.5 16])
title('QE scheme for small \nu')
%xlabel('log(nu)')
xlabel('\nu')
ylabel('Option price')
hold on
plot(nu(1:n),v0,'--g','LineWidth',1)
% hold on
% plot(nu(1:n),QEzero,'--g','LineWidth',2)
hold off
legend('QE','Modified QE', 'B-S')
legend('location','southeast')
%}
%% Generate LaTex code for a table of "nu is close or equal to zero"
clear input;
temp = NaN(1,8);
call_exact = [ourGBMCallPrice.exactPrice,temp];%,SobolCallPrice_exact,temp];
%put_exact = [GBMPutPrice_exact,temp];%,SobolPutPrice_exact,temp];
content = [call_exact; QEprice; ourQEprice;dif];
% content = [call_exact;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_exact;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'0','1e-6','1e-5','1e-4','1e-3','1e-2','0.05','0.1','0.5'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'Exact price','QE','Modified QE','Diff'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat={'%.3f',9};%,'%.4f',3,'%5.2fs',3,'%3.2E',3};
input.dataNanString = '-';
input.tableColumnAlignment = 'c';
% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = '$\nu$ is close or equal to zero';
% LaTex table label:
input.tableLabel = 'nu=0';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Switch to landscape table:
input.landscape = 0;
% call latexTable:
latex = latexTable(input);