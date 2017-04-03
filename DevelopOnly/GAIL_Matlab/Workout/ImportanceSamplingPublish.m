%{
%% Importance Sampling
%
% Set assetPath parameters
clearvars
T=5;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector
inp.assetParam.interest = 0.0;          % risk-free interest rate
inp.assetParam.volatility = 0.3;        % fixed vlatility of asset prices
inp.assetParam.Vinst = 0.5;%0.09            % initial value of volatility
inp.assetParam.Vlong = 0.16;%0.09            % theta
inp.assetParam.kappa = 1;               % kappa
inp.assetParam.nu = 0.4;%0                  % volatility of asset price volatility
inp.assetParam.rho = -0.3;%0                 % rho
inp.assetParam.pathType = 'GBM'; 
inp.priceParam.cubMethod = 'IID_MC';
%Set optPayoff parameter
inp.assetParam.initPrice = 60;            % strike price
%Set error tolerance
inp.priceParam.absTol = 0;              % absolute tolerance
inp.priceParam.relTol = 0.01;           % one penny on the dollar relative tolerance
strike = [20,60,100];

%% European call option
%
% * cubMethod = 'IID_MC'
%With importance sampling
GBMCallPrice_withIS = zeros(1,3);
GBMCallPriceIS_paths = zeros(1,3);
GBMCallPriceIS_time = zeros(1,3);
GBMCallPrice_exact = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike >= inp.assetParam.initPrice
        inp.assetParam.meanShift = 5;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    ourGBMCallPrice = optPrice(inp);
    [GBMCallPrice_withIS(i), out] = genOptPrice(ourGBMCallPrice);
    GBMCallPriceIS_paths(i) = out.nPaths;
    GBMCallPriceIS_time(i) = out.time;
    GBMCallPrice_exact(i) = ourGBMCallPrice.exactPrice;
end
iid_temp = [GBMCallPrice_withIS,GBMCallPriceIS_time,GBMCallPriceIS_paths];
% GBMCP_exact = num2str(GBMCallPrice_exact(:),'$%.4f \t')
GBMCPI = sprintf('$%.4f \n',GBMCallPrice_withIS);
GBMCP_exact = sprintf('$%.4f \t',GBMCallPrice_exact); 
GBMCPI_paths={};
GBMCPI_paths =[GBMCPI_paths {sprintf('%3.2E\n',GBMCallPriceIS_paths(:))}];
GBMCPI_time = sprintf('%5.2f\t',GBMCallPriceIS_time);
%strike;
GBMCP_exact;
GBMCPI;
GBMCPI_paths
% strrep(GBMCPI_paths,'e+00','e+')
GBMCPI_time;
%return
% Without importance Sampling
inp.assetParam.meanShift = 0;
GBMCallPrice_withoutIS=zeros(1,3);
GBMCallPrice_paths = zeros(1,3);
GBMCallPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourGBMCallPrice = optPrice(inp);
    [GBMCallPrice_withoutIS(i), out] = genOptPrice(ourGBMCallPrice);
    GBMCallPrice_paths(i) = out.nPaths;
    GBMCallPrice_time(i) = out.time;
end
iid_temp = [iid_temp;GBMCallPrice_withoutIS,GBMCallPrice_time,GBMCallPrice_paths];
GBMCP = sprintf('$%.4f \t',GBMCallPrice_withoutIS);
GBMCP_paths = sprintf('%3.2E \t',GBMCallPrice_paths);
GBMCP_time = sprintf('%5.2fs\t',GBMCallPrice_time);
strike;
GBMCP;
GBMCP_paths;
GBMCP_time;
%%
% * cubMethod = 'Sobol'
%Set assetPath parameters
T=5;
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector       
inp.priceParam.cubMethod = 'Sobol';

%With importance sampling
SobolCallPrice_withIS=zeros(1,3);
SobolCallPriceIS_paths = zeros(1,3);
SobolCallPriceIS_time = zeros(1,3);
SobolCallPrice_exact = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike >= inp.assetParam.initPrice
        inp.assetParam.meanShift = 5;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    ourSobolCallPrice = optPrice(inp);
    [SobolCallPrice_withIS(i), out] = genOptPrice(ourSobolCallPrice);
    SobolCallPriceIS_paths(i) = out.nPaths;
    SobolCallPriceIS_time(i) = out.time;
    SobolCallPrice_exact(i) = ourSobolCallPrice.exactPrice;
end
Sobol_temp = [SobolCallPrice_withIS,SobolCallPriceIS_time,SobolCallPriceIS_paths];
SobolCPI = sprintf('$%.4f \t',SobolCallPrice_withIS);
SobolCP_exact = sprintf('$%.4f \t',SobolCallPrice_exact); 
SobolCPI_paths = sprintf('%3.2E \t',SobolCallPriceIS_paths);
SobolCPI_time = sprintf('%5.2fs\t',SobolCallPriceIS_time);
strike;
SobolCP_exact;
SobolCPI;
SobolCPI_paths;
SobolCPI_time;

% Without importance Sampling
inp.assetParam.meanShift = 0;
SobolCallPrice_withoutIS=zeros(1,3);
SobolCallPrice_paths = zeros(1,3);
SobolCallPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourSobolCallPrice = optPrice(inp);
    [SobolCallPrice_withoutIS(i), out] = genOptPrice(ourSobolCallPrice);
    SobolCallPrice_paths(i) = out.nPaths;
    SobolCallPrice_time(i) = out.time;
end
Sobol_temp = [Sobol_temp;SobolCallPrice_withoutIS,SobolCallPrice_time,SobolCallPrice_paths];
SobolCP = sprintf('$%.4f \t',SobolCallPrice_withoutIS);
SobolCP_paths = sprintf('%3.2E \t',SobolCallPrice_paths);
SobolCP_time = sprintf('%5.2fs\t',SobolCallPrice_time);
strike;
SobolCP;
SobolCP_paths;
SobolCP_time;

%% European put option

T=5;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector
inp.assetParam.interest = 0.0;          % risk-free interest rate
inp.assetParam.volatility = 0.3;        % fixed vlatility of asset prices
inp.assetParam.Vinst = 0.09;            % initial value of volatility
inp.assetParam.Vlong = 0.09;            % theta
inp.assetParam.kappa = 1;               % kappa
inp.assetParam.nu = 0;                  % volatility of asset price volatility
inp.assetParam.rho = 0;                 % rho
inp.assetParam.pathType = 'GBM'; 
inp.priceParam.cubMethod = 'IID_MC';

inp.payoffParam.putCallType = {'put'};
%strike = fliplr(strike);   %initialPrice = [40,60,80]
%%
% * cubMethod = 'IID_MC'
%With Importance Sampling
GBMPutPrice_withIS=zeros(1,3);
GBMPutPriceIS_paths = zeros(1,3);
GBMPutPriceIS_time = zeros(1,3);
GBMPutPrice_exact = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike <= inp.assetParam.initPrice
        inp.assetParam.meanShift = -2;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    ourGBMPutPrice = optPrice(inp);
    [GBMPutPrice_withIS(i), out] = genOptPrice(ourGBMPutPrice);
    GBMPutPriceIS_paths(i) = out.nPaths;
    GBMPutPriceIS_time(i) = out.time;
    GBMPutPrice_exact(i) = ourGBMPutPrice.exactPrice;
end
iid_temp = [iid_temp;GBMPutPrice_withIS,GBMPutPriceIS_time,GBMPutPriceIS_paths];
GBMPPI = sprintf('$%.4f \t',GBMPutPrice_withIS);
GBMPP_exact = sprintf('$%.4f \t',GBMPutPrice_exact); 
GBMPPI_paths = sprintf('%3.2E \t',GBMPutPriceIS_paths);
GBMPPI_time = sprintf('%5.2fs\t',GBMPutPriceIS_time);
strike;
GBMPP_exact;
GBMPPI;
GBMPPI_paths;
GBMPPI_time;

% Without Importance Sampling
inp.assetParam.meanShift = 0;
GBMPutPrice_withoutIS=zeros(1,3);
GBMPutPrice_paths = zeros(1,3);
GBMPutPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourGBMPutPrice = optPrice(inp);
    [GBMPutPrice_withoutIS(i), out] = genOptPrice(ourGBMPutPrice);
    GBMPutPrice_paths(i) = out.nPaths;
    GBMPutPrice_time(i) = out.time;
end
iid_temp = [iid_temp;GBMPutPrice_withoutIS,GBMPutPrice_time,GBMPutPrice_paths];
GBMPP = sprintf('$%.4f \t',GBMPutPrice_withoutIS);
GBMPP_paths = sprintf('%3.2E \t',GBMPutPrice_paths);
GBMPP_time = sprintf('%5.2fs\t',GBMPutPrice_time);
strike;
GBMPP;
GBMPP_paths;
GBMPP_time;
%%
% * cubMethod = 'Sobol'
%Set assetPath parameters
T=5;
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector       
inp.priceParam.cubMethod = 'Sobol';

%With importance sampling
SobolPutPrice_withIS=zeros(1,3);
SobolPutPriceIS_paths = zeros(1,3);
SobolPutPriceIS_time = zeros(1,3);
SobolPutPrice_exact = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike <= inp.assetParam.initPrice
        inp.assetParam.meanShift = -2;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    ourSobolPutPrice = optPrice(inp);
    [SobolPutPrice_withIS(i), out] = genOptPrice(ourSobolPutPrice);
    SobolPutPriceIS_paths(i) = out.nPaths;
    SobolPutPriceIS_time(i) = out.time;
    SobolPutPrice_exact(i) = ourSobolPutPrice.exactPrice;
end
Sobol_temp = [Sobol_temp;SobolPutPrice_withIS,SobolPutPriceIS_time,SobolPutPriceIS_paths];
SobolPPI = sprintf('$%.4f \t',SobolPutPrice_withIS);
SobolPP_exact = sprintf('$%.4f \t',SobolPutPrice_exact); 
SobolPPI_paths = sprintf('%3.2E \t',SobolPutPriceIS_paths);
SobolPPI_time = sprintf('%5.2fs\t',SobolPutPriceIS_time);
strike;
SobolPP_exact;
SobolPPI;
SobolPPI_paths;
SobolPPI_time;

% Without importance Sampling
inp.assetParam.meanShift = 0;
SobolPutPrice_withoutIS=zeros(1,3);
SobolPutPrice_paths = zeros(1,3);
SobolPutPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourSobolPutPrice = optPrice(inp);
    [SobolPutPrice_withoutIS(i), out] = genOptPrice(ourSobolPutPrice);
    SobolPutPrice_paths(i) = out.nPaths;
    SobolPutPrice_time(i) = out.time;
end
Sobol_temp = [Sobol_temp;SobolPutPrice_withoutIS,SobolPutPrice_time,SobolPutPrice_paths];
SobolPP = sprintf('$%.4f \t',SobolPutPrice_withoutIS);
SobolPP_paths = sprintf('%3.2E \t',SobolPutPrice_paths);
SobolPP_time = sprintf('%5.2fs\t',SobolPutPrice_time);
strike;
SobolPP;
SobolPP_paths;
SobolPP_time;

%% American put option
% * cubMethod = 'IID_MC'
% American put option with importance sampling
inp.payoffParam.optType = {'american'}; %change from European to American
inp.priceParam.cubMethod = 'IID_MC';
AmericanPutPrice_withIS=zeros(1,3);
AmericanPutPriceIS_paths = zeros(1,3);
AmericanPutPriceIS_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike <= inp.assetParam.initPrice
        inp.assetParam.meanShift = -2;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    [AmericanPutPrice_withIS(i), out] = genOptPrice(AmericanPut);
    AmericanPutPriceIS_paths(i) = out.nPaths;
    AmericanPutPriceIS_time(i) = out.time;   
end
APPI = sprintf('$%.4f \t',AmericanPutPrice_withIS);
APPI_paths = sprintf('%3.2E \t',AmericanPutPriceIS_paths);
APPI_time = sprintf('%5.2fs\t',AmericanPutPriceIS_time);
strike;
APPI;
APPI_paths;
APPI_time;
iid_temp = [iid_temp;AmericanPutPrice_withIS,AmericanPutPriceIS_time,...
    AmericanPutPriceIS_paths];
% American put option without importance sampling
inp.assetParam.meanShift = 0;
AmericanPutPrice_withoutIS=zeros(1,3);
AmericanPutPrice_paths = zeros(1,3);
AmericanPutPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    [AmericanPutPrice_withoutIS(i), out] = genOptPrice(AmericanPut);
    AmericanPutPrice_paths(i) = out.nPaths;
    AmericanPutPrice_time(i) = out.time; 
end
APP = sprintf('$%.4f \t',AmericanPutPrice_withoutIS);
APP_paths = sprintf('%3.2E \t',AmericanPutPrice_paths);
APP_time = sprintf('%5.2fs\t',AmericanPutPrice_time);
strike;
APP;
APP_paths;
APP_time;
iid_temp = [iid_temp;AmericanPutPrice_withoutIS,AmericanPutPrice_time,...
    AmericanPutPrice_paths];
%%
% * cubMethod = 'Sobol'
%Set assetPath parameters
T=5;
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector       
inp.priceParam.cubMethod = 'Sobol';

%With importance sampling
AmericanPutPrice_withIS=zeros(1,3);
AmericanPutPriceIS_paths = zeros(1,3);
AmericanPutPriceIS_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike <= inp.assetParam.initPrice
        inp.assetParam.meanShift = -2;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    [AmericanPutPrice_withIS(i), out] = genOptPrice(AmericanPut);
    AmericanPutPriceIS_paths(i) = out.nPaths;
    AmericanPutPriceIS_time(i) = out.time;   
end
Sobol_temp = [Sobol_temp;AmericanPutPrice_withIS,AmericanPutPriceIS_time,...
    AmericanPutPriceIS_paths];
%Without importance sampling
inp.assetParam.meanShift = 0;
AmericanPutPrice_withIS=zeros(1,3);
AmericanPutPriceIS_paths = zeros(1,3);
AmericanPutPriceIS_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    [AmericanPutPrice_withIS(i), out] = genOptPrice(AmericanPut);
    AmericanPutPriceIS_paths(i) = out.nPaths;
    AmericanPutPriceIS_time(i) = out.time;   
end
Sobol_temp = [Sobol_temp;AmericanPutPrice_withIS,AmericanPutPriceIS_time,...
    AmericanPutPriceIS_paths];

%}
%% Generate LaTex code for a table of "call and put options",pathtype='GBM',cubmethod='IID_MC'
clear input;
temp = NaN(1,6);
call_exact = [GBMCallPrice_exact,temp];%,SobolCallPrice_exact,temp];
put_exact = [GBMPutPrice_exact,temp];%,SobolPutPrice_exact,temp];
content = [call_exact;iid_temp(1:2,:); put_exact;iid_temp(3:end,:)];
% content = [call_exact;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_exact;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'K=20','K=60','K=100','K=20','K=60','K=100','K=20','K=60','K=100'};
% Set row labels (use empty string for no label):
% input.tableRowLabels = {'EuropeanCall\_exact price','EuropeanCall\_withIS',...
% 'EuropeanCall\_withoutIS','EuropeanPut\_exact price','EuropeanPut\_withIS',...
% 'EuropeanPut\_withoutIS','AmericanPut\_withIS','AmericanPut\_withoutIS'};
input.tableRowLabels = {'ECexact','ECwIS',...
'ECw/oIS','EPexact','EPwIS',...
'EPw/oIS','APwIS','APw/oIS'};
% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat={'%.3f',3,'%5.2f',3,'%3.2e',3};%,'%.4f',3,'%5.2fs',3,'%3.2E',3};
input.dataNanString = '-';
input.tableColumnAlignment = 'c';
% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 0;
% LaTex table caption:
input.tableCaption = 'IID\_MC S0=60';
% LaTex table label:
input.tableLabel = 'MyTableLabel';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Switch to landscape table:
input.landscape = 0;
% call latexTable:
latex = latexTable(input);

%% Generate LaTex code for a table of "call and put options",pathtype='GBM',cubmethod='Sobol'
clear input;
temp = NaN(1,6);
call_exact = [GBMCallPrice_exact,temp];%,SobolCallPrice_exact,temp];
put_exact = [GBMPutPrice_exact,temp];%,SobolPutPrice_exact,temp];
content = [call_exact;Sobol_temp(1:2,:); put_exact;Sobol_temp(3:end,:)];
% content = [call_exact;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_exact;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'K=20','K=60','K=100','K=20','K=60','K=100','K=20','K=60','K=100'};
% Set row labels (use empty string for no label):
% input.tableRowLabels = {'EuropeanCall\_exact price','EuropeanCall\_withIS',...
% 'EuropeanCall\_withoutIS','EuropeanPut\_exact price','EuropeanPut\_withIS',...
% 'EuropeanPut\_withoutIS','AmericanPut\_withIS','AmericanPut\_withoutIS'};
input.tableRowLabels = {'ECexact','ECwIS',...
'ECw/oIS','EPexact','EPwIS',...
'EPw/oIS','APwIS','APw/oIS'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used
input.dataFormat={'%.3f',3,'%5.2f',3,'%3.2e',3};%,'%.4f',3,'%5.2fs',3,'%3.2E',3};
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

%% Reference

