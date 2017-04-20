%%{
%% Importance Sampling
%
% Set assetPath parameters
clearvars
T=5;
delta_t=0.2;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector
inp.assetParam.interest = 0.0;          % risk-free interest rate
inp.assetParam.volatility = 0.4;        % fixed vlatility of asset prices
inp.assetParam.Vinst = 0.5;            % initial value of volatility
inp.assetParam.Vlong = 0.16;            % theta
inp.assetParam.kappa = 1;               % kappa
inp.assetParam.nu = 0.4;                % volatility of asset price volatility
inp.assetParam.rho = -0.3;              % rho
inp.assetParam.pathType = 'QE'; 
inp.priceParam.cubMethod = 'IID_MC';
%Set optPayoff parameter
inp.assetParam.initPrice = 60;            % strike price
%Set error tolerance
inp.priceParam.absTol = 0;              % absolute tolerance
inp.priceParam.relTol = 0.01;           % one penny on the dollar relative tolerance
strike = [20,60,100];
avg = 5; % run five times and take average
%{
%% * cubMethod = 'IID_MC'
%*****************************
% European call option
%*****************************

%With importance sampling
IIDCallPrice_withIS = zeros(1,3);
IIDCallPriceIS_paths = zeros(1,3);
IIDCallPriceIS_time = zeros(1,3);
IIDCallPrice_QE = zeros(1,3);
IIDCallElapsed = zeros(1,3);
inp.assetParam.meanShift = 2;
for i=1:3
    inp.payoffParam.strike = strike(i);
%     if inp.payoffParam.strike >= inp.assetParam.initPrice
%         inp.assetParam.meanShift = 2;
%     else
%         inp.assetParam.meanShift = 0;
%     end
    %Construct an optPrice object
    ourIIDCallPrice = optPrice(inp);
    for j = 1:avg        
        [temp, out] = genOptPrice(ourIIDCallPrice);
        IIDCallPrice_withIS(i) = IIDCallPrice_withIS(i) + temp;
        IIDCallPriceIS_paths(i) = IIDCallPriceIS_paths(i) + out.nPaths;
        IIDCallPriceIS_time(i) = IIDCallPriceIS_time(i) + out.time;
        
        % Calculate option price by provided codes
        %  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
        tic
        Ntime = T/delta_t; 
        [a,b] = MC_QE(inp.assetParam.initPrice,inp.assetParam.interest,0,T,...
            inp.assetParam.Vinst,inp.assetParam.Vlong,inp.assetParam.kappa,inp.assetParam.nu,...
            inp.assetParam.rho,Ntime,1e6,1);
        PT = a(:,Ntime + 1);
        PT = max(PT-inp.payoffParam.strike,0);
        PP = mean(PT);
        temp = PP*exp(-inp.assetParam.interest*T);
        IIDCallPrice_QE(i) = IIDCallPrice_QE(i) + temp;
        IIDCallElapsed(i) = IIDCallElapsed(i) + toc;
    end
end
IIDCallPrice_withIS = IIDCallPrice_withIS ./ avg;
IIDCallPriceIS_paths = IIDCallPriceIS_paths ./ avg;
IIDCallPriceIS_time = IIDCallPriceIS_time ./ avg;

IIDCallPrice_QE = IIDCallPrice_QE ./ avg;
IIDCallElapsed = IIDCallElapsed ./ avg;

iid_temp = [IIDCallPrice_withIS,IIDCallPriceIS_time,IIDCallPriceIS_paths];
%IIDCP_QE = num2str(IIDCallPrice_QE(:),'$%.4f \t');
IIDCPI = sprintf('$%.4f \n',IIDCallPrice_withIS);
IIDCP_QE = sprintf('$%.4f \t',IIDCallPrice_QE); 
IIDCPI_paths = sprintf('%3.2E \t',IIDCallPriceIS_paths);
IIDCPI_time = sprintf('%5.2fs\t',IIDCallPriceIS_time);
strike;
IIDCP_QE;
IIDCPI;
IIDCPI_paths;
IIDCPI_time;

% Without importance Sampling
inp.assetParam.meanShift = 0;
IIDCallPrice_withoutIS=zeros(1,3);
IIDCallPrice_paths = zeros(1,3);
IIDCallPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourIIDCallPrice = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(ourIIDCallPrice);
        IIDCallPrice_withoutIS(i) = IIDCallPrice_withoutIS(i) + temp;
        IIDCallPrice_paths(i) = IIDCallPrice_paths(i) + out.nPaths;
        IIDCallPrice_time(i) = IIDCallPrice_time(i) + out.time;
    end
end
IIDCallPrice_withoutIS = IIDCallPrice_withoutIS ./ avg;
IIDCallPrice_paths = IIDCallPrice_paths ./ avg;
IIDCallPrice_time = IIDCallPrice_time ./ avg;

iid_temp = [iid_temp;IIDCallPrice_withoutIS,IIDCallPrice_time,IIDCallPrice_paths];
IIDCP = sprintf('$%.4f \t',IIDCallPrice_withoutIS);
IIDCP_paths = sprintf('%3.2E \t',IIDCallPrice_paths);
IIDCP_time = sprintf('%5.2fs\t',IIDCallPrice_time);
strike;
IIDCP;
IIDCP_paths;
IIDCP_time;

%************************
% European put option
%************************
% T=5;
% delta_t=0.2;
% t0 = delta_t;
% inp.timeDim.timeVector = t0:delta_t:T;  % time vector
% inp.assetParam.interest = 0.0;          % risk-free interest rate
% inp.assetParam.volatility = 0.3;        % fixed vlatility of asset prices
% inp.assetParam.Vinst = 0.09;            % initial value of volatility
% inp.assetParam.Vlong = 0.09;            % theta
% inp.assetParam.kappa = 1;               % kappa
% inp.assetParam.nu = 0;                  % volatility of asset price volatility
% inp.assetParam.rho = 0;                 % rho
% inp.assetParam.pathType = 'GBM'; 
% inp.priceParam.cubMethod = 'IID_MC';
inp.payoffParam.putCallType = {'put'};
% inp.priceParam.relTol = 0.01;    
%strike = fliplr(strike);   %initialPrice = [40,60,80]
%
%
%With Importance Sampling
IIDPutPrice_withIS=zeros(1,3);
IIDPutPriceIS_paths = zeros(1,3);
IIDPutPriceIS_time = zeros(1,3);
IIDPutPrice_QE = zeros(1,3);
IIDPutElapsed = zeros(1,3);
inp.assetParam.meanShift = -0.8;
for i=1:3
    inp.payoffParam.strike = strike(i);
%     if inp.payoffParam.strike <= inp.assetParam.initPrice
%         inp.assetParam.meanShift = -0.8;
%     else
%         inp.assetParam.meanShift = 0;
%     end
    %Construct an optPrice object
    ourIIDPutPrice = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(ourIIDPutPrice);
        IIDPutPrice_withIS(i) = IIDPutPrice_withIS(i) + temp;
        IIDPutPriceIS_paths(i) = IIDPutPriceIS_paths(i) + out.nPaths;
        IIDPutPriceIS_time(i) = IIDPutPriceIS_time(i) + out.time;
        
        % Calculate option price by provided codes
        %  MC_QE(S0,r,d,T,Vinst,Vlong,kappa,epsilon,rho,NTime,NSim,NBatches)
        tic
        Ntime = T/delta_t; 
        [a,b] = MC_QE(inp.assetParam.initPrice,inp.assetParam.interest,0,T,...
            inp.assetParam.Vinst,inp.assetParam.Vlong,inp.assetParam.kappa,inp.assetParam.nu,...
            inp.assetParam.rho,Ntime,1e6,1);
        PT = a(:,Ntime + 1);
        PT = max(inp.payoffParam.strike-PT,0);
        PP = mean(PT);
        temp = PP*exp(-inp.assetParam.interest*T);
        IIDPutPrice_QE(i) = IIDPutPrice_QE(i) + temp;
        IIDPutElapsed(i) = IIDPutElapsed(i) + toc;
    end
end
IIDPutPrice_withIS = IIDPutPrice_withIS ./ avg;
IIDPutPriceIS_paths = IIDPutPriceIS_paths ./ avg;
IIDPutPriceIS_time = IIDPutPriceIS_time ./ avg;

IIDPutPrice_QE = IIDPutPrice_QE ./ avg;
IIDPutElapsed = IIDPutElapsed ./ avg;

iid_temp = [iid_temp;IIDPutPrice_withIS,IIDPutPriceIS_time,IIDPutPriceIS_paths];
IIDPPI = sprintf('$%.4f \t',IIDPutPrice_withIS);
IIDPP_QE = sprintf('$%.4f \t',IIDPutPrice_QE); 
IIDPPI_paths = sprintf('%3.2E \t',IIDPutPriceIS_paths);
IIDPPI_time = sprintf('%5.2fs\t',IIDPutPriceIS_time);
strike;
IIDPP_QE;
IIDPPI;
IIDPPI_paths;
IIDPPI_time;

% Without Importance Sampling
inp.assetParam.meanShift = 0;
IIDPutPrice_withoutIS=zeros(1,3);
IIDPutPrice_paths = zeros(1,3);
IIDPutPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourIIDPutPrice = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(ourIIDPutPrice);
        IIDPutPrice_withoutIS(i) = IIDPutPrice_withoutIS(i) + temp;
        IIDPutPrice_paths(i) = IIDPutPrice_paths(i) + out.nPaths;
        IIDPutPrice_time(i) = IIDPutPrice_time(i) + out.time;
    end
end
IIDPutPrice_withoutIS = IIDPutPrice_withoutIS ./ avg;
IIDPutPrice_paths = IIDPutPrice_paths ./ avg;
IIDPutPrice_time = IIDPutPrice_time ./ avg;

iid_temp = [iid_temp;IIDPutPrice_withoutIS,IIDPutPrice_time,IIDPutPrice_paths];
IIDPP = sprintf('$%.4f \t',IIDPutPrice_withoutIS);
IIDPP_paths = sprintf('%3.2E \t',IIDPutPrice_paths);
IIDPP_time = sprintf('%5.2fs\t',IIDPutPrice_time);
strike;
IIDPP;
IIDPP_paths;
IIDPP_time;
%}
%*****************************
% American put option
%*****************************
% American put option with importance sampling
% delta_t=0.2;
% t0 = delta_t;
% inp.timeDim.timeVector = t0:delta_t:T;
inp.payoffParam.optType = {'american'}; %change from European to American
inp.payoffParam.putCallType = {'put'};
% inp.priceParam.cubMethod = 'IID_MC';
% inp.priceParam.relTol = 0.01;    

AmericanPutPrice_withIS=zeros(1,3);
AmericanPutPriceIS_paths = zeros(1,3);
AmericanPutPriceIS_time = zeros(1,3);
% inp.assetParam.meanShift = -0.8;
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike < inp.assetParam.initPrice
        inp.assetParam.meanShift = -0.8;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    for j = 1: avg        
        [temp, out] = genOptPrice(AmericanPut);
        AmericanPutPrice_withIS(i) = AmericanPutPrice_withIS(i) + temp;
        AmericanPutPriceIS_paths(i) = AmericanPutPriceIS_paths(i) + out.nPaths;
        AmericanPutPriceIS_time(i) = AmericanPutPriceIS_time(i) + out.time;   
    end
end
AmericanPutPrice_withIS = AmericanPutPrice_withIS ./ avg;
AmericanPutPriceIS_paths = AmericanPutPriceIS_paths ./ avg;
AmericanPutPriceIS_time = AmericanPutPriceIS_time ./ avg;

APPI = sprintf('$%.4f \t',AmericanPutPrice_withIS);
APPI_paths = sprintf('%3.2E \t',AmericanPutPriceIS_paths);
APPI_time = sprintf('%5.2fs\t',AmericanPutPriceIS_time);
strike;
APPI;
APPI_paths;
APPI_time;
% iid_temp = [iid_temp;AmericanPutPrice_withIS,AmericanPutPriceIS_time,...
%     AmericanPutPriceIS_paths];

% American put option without importance sampling
inp.assetParam.meanShift = 0;
AmericanPutPrice_withoutIS=zeros(1,3);
AmericanPutPrice_paths = zeros(1,3);
AmericanPutPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(AmericanPut);
        AmericanPutPrice_withoutIS(i) = AmericanPutPrice_withoutIS(i) + temp;
        AmericanPutPrice_paths(i) = AmericanPutPrice_paths(i) + out.nPaths;
        AmericanPutPrice_time(i) = AmericanPutPrice_time(i) + out.time; 
    end
end
AmericanPutPrice_withoutIS = AmericanPutPrice_withoutIS ./ avg;
AmericanPutPrice_paths = AmericanPutPrice_paths ./ avg;
AmericanPutPrice_time = AmericanPutPrice_time ./ avg;

APP = sprintf('$%.4f \t',AmericanPutPrice_withoutIS);
APP_paths = sprintf('%3.2E \t',AmericanPutPrice_paths);
APP_time = sprintf('%5.2fs\t',AmericanPutPrice_time);
strike;
APP;
APP_paths;
APP_time;
% iid_temp = [iid_temp;AmericanPutPrice_withoutIS,AmericanPutPrice_time,...
%     AmericanPutPrice_paths];
%}
return
%{
%% * cubMethod = 'Sobol'

%*****************************
% European call option
%*****************************
%Set assetPath parameters
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector 
inp.payoffParam.optType = {'euro'};
inp.payoffParam.putCallType = {'call'};
inp.priceParam.cubMethod = 'Sobol';
inp.priceParam.relTol = 0.001;    
%With importance sampling
SobolCallPrice_withIS=zeros(1,3);
SobolCallPriceIS_paths = zeros(1,3);
SobolCallPriceIS_time = zeros(1,3);
%SobolCallPrice_QE = zeros(1,3);
% inp.assetParam.meanShift = 1;
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike >= inp.assetParam.initPrice
        inp.assetParam.meanShift = 1;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    ourSobolCallPrice = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(ourSobolCallPrice);
        SobolCallPrice_withIS(i) = SobolCallPrice_withIS(i) + temp;
        SobolCallPriceIS_paths(i) = SobolCallPriceIS_paths(i) + out.nPaths;
        SobolCallPriceIS_time(i) = SobolCallPriceIS_time(i) + out.time;
        %SobolCallPrice_QE(i) = ourSobolCallPrice.QEPrice;
    end
end
SobolCallPrice_withIS = SobolCallPrice_withIS ./ avg;
SobolCallPriceIS_paths = SobolCallPriceIS_paths ./ avg;
SobolCallPriceIS_time = SobolCallPriceIS_time ./ avg;
    
Sobol_temp = [SobolCallPrice_withIS,SobolCallPriceIS_time,SobolCallPriceIS_paths];
SobolCPI = sprintf('$%.4f \t',SobolCallPrice_withIS);
%SobolCP_QE = sprintf('$%.4f \t',SobolCallPrice_QE); 
SobolCPI_paths = sprintf('%3.2E \t',SobolCallPriceIS_paths);
SobolCPI_time = sprintf('%5.2fs\t',SobolCallPriceIS_time);
strike;
%SobolCP_QE;
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
    for j = 1:avg
        [temp, out] = genOptPrice(ourSobolCallPrice);
        SobolCallPrice_withoutIS(i) = SobolCallPrice_withoutIS(i) + temp;
        SobolCallPrice_paths(i) = SobolCallPrice_paths(i) + out.nPaths;
        SobolCallPrice_time(i) = SobolCallPrice_time(i) + out.time;
    end 
end
SobolCallPrice_withoutIS = SobolCallPrice_withoutIS ./ avg;
SobolCallPrice_paths = SobolCallPrice_paths ./avg;
SobolCallPrice_time = SobolCallPrice_time ./ avg;

Sobol_temp = [Sobol_temp;SobolCallPrice_withoutIS,SobolCallPrice_time,SobolCallPrice_paths];
SobolCP = sprintf('$%.4f \t',SobolCallPrice_withoutIS);
SobolCP_paths = sprintf('%3.2E \t',SobolCallPrice_paths);
SobolCP_time = sprintf('%5.2fs\t',SobolCallPrice_time);
strike;
SobolCP;
SobolCP_paths;
SobolCP_time;
%%}

%************************
% European put option
%************************
%Set assetPath parameters
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector       
% inp.priceParam.cubMethod = 'Sobol';
inp.payoffParam.putCallType = {'put'};
% inp.priceParam.relTol = 0.001;    
%With importance sampling
SobolPutPrice_withIS=zeros(1,3);
SobolPutPriceIS_paths = zeros(1,3);
SobolPutPriceIS_time = zeros(1,3);
%inp.assetParam.meanShift = -1;
for i=1:3
    inp.payoffParam.strike = strike(i);
    if inp.payoffParam.strike < inp.assetParam.initPrice
        inp.assetParam.meanShift = -1;
    else
        inp.assetParam.meanShift = 0;
    end
    %Construct an optPrice object
    ourSobolPutPrice = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(ourSobolPutPrice);
        SobolPutPrice_withIS(i) = SobolPutPrice_withIS(i) + temp;
        SobolPutPriceIS_paths(i) = SobolPutPriceIS_paths(i) + out.nPaths;
        SobolPutPriceIS_time(i) = SobolPutPriceIS_time(i) + out.time;
    end
end
SobolPutPrice_withIS = SobolPutPrice_withIS ./ avg;
SobolPutPriceIS_paths = SobolPutPriceIS_paths ./ avg;
SobolPutPriceIS_time = SobolPutPriceIS_time ./ avg;

Sobol_temp = [Sobol_temp;SobolPutPrice_withIS,SobolPutPriceIS_time,SobolPutPriceIS_paths];
SobolPPI = sprintf('$%.4f \t',SobolPutPrice_withIS);
%SobolPP_QE = sprintf('$%.4f \t',SobolPutPrice_QE); 
SobolPPI_paths = sprintf('%3.2E \t',SobolPutPriceIS_paths);
SobolPPI_time = sprintf('%5.2fs\t',SobolPutPriceIS_time);
strike;
%SobolPP_QE;
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
    for j = 1: avg
        [temp, out] = genOptPrice(ourSobolPutPrice);
        SobolPutPrice_withoutIS(i) = SobolPutPrice_withoutIS(i) + temp;
        SobolPutPrice_paths(i) = SobolPutPrice_paths(i) + out.nPaths;
        SobolPutPrice_time(i) = SobolPutPrice_time(i) + out.time;
    end
end
SobolPutPrice_withoutIS = SobolPutPrice_withoutIS ./ avg;
SobolPutPrice_paths = SobolPutPrice_paths ./ avg; 
SobolPutPrice_time = SobolPutPrice_time ./ avg; 

Sobol_temp = [Sobol_temp;SobolPutPrice_withoutIS,SobolPutPrice_time,SobolPutPrice_paths];
SobolPP = sprintf('$%.4f \t',SobolPutPrice_withoutIS);
SobolPP_paths = sprintf('%3.2E \t',SobolPutPrice_paths);
SobolPP_time = sprintf('%5.2fs\t',SobolPutPrice_time);
strike;
SobolPP;
SobolPP_paths;
SobolPP_time;



%*****************************
% American put option
%*****************************
%Set assetPath parameters
delta_t=0.5;
t0 = delta_t;
inp.timeDim.timeVector = t0:delta_t:T;  % time vector       
inp.priceParam.cubMethod = 'Sobol';
inp.payoffParam.optType = {'american'};
% inp.payoffParam.putCallType = {'put'};
% inp.priceParam.relTol = 0.001;    
%With importance sampling
AmericanPutPrice_withIS=zeros(1,3);
AmericanPutPriceIS_paths = zeros(1,3);
AmericanPutPriceIS_time = zeros(1,3);
inp.assetParam.meanShift = -0.3;
for i=1:3
    inp.payoffParam.strike = strike(i);
%     if inp.payoffParam.strike <= inp.assetParam.initPrice
%         inp.assetParam.meanShift = -0.8;
%     else
%         inp.assetParam.meanShift = 0;
%     end
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    for j = 1: avg
        [temp, out] = genOptPrice(AmericanPut);
        AmericanPutPrice_withIS(i) = AmericanPutPrice_withIS(i) + temp;
        AmericanPutPriceIS_paths(i) = AmericanPutPriceIS_paths(i) + out.nPaths;
        AmericanPutPriceIS_time(i) = AmericanPutPriceIS_time(i) + out.time;   
    end
end
 AmericanPutPrice_withIS = AmericanPutPrice_withIS ./ avg;
 AmericanPutPriceIS_paths = AmericanPutPriceIS_paths ./ avg;
 AmericanPutPriceIS_time = AmericanPutPriceIS_time ./ avg;

Sobol_temp = [Sobol_temp;AmericanPutPrice_withIS,AmericanPutPriceIS_time,...
    AmericanPutPriceIS_paths];
%Without importance sampling
inp.assetParam.meanShift = 0;
AmericanPutPrice_withoutIS=zeros(1,3);
AmericanPutPrice_paths = zeros(1,3);
AmericanPutPrice_time = zeros(1,3);
for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    AmericanPut = optPrice(inp);
    for j = 1:avg
        [temp, out] = genOptPrice(AmericanPut);
        AmericanPutPrice_withoutIS(i) = AmericanPutPrice_withoutIS(i) + temp;
        AmericanPutPrice_paths(i) = AmericanPutPrice_paths(i) + out.nPaths;
        AmericanPutPrice_time(i) = AmericanPutPrice_time(i) + out.time;   
    end
end
AmericanPutPrice_withoutIS = AmericanPutPrice_withoutIS ./ avg;
AmericanPutPrice_paths = AmericanPutPrice_paths ./ avg;
AmericanPutPrice_time = AmericanPutPrice_time ./ avg;

Sobol_temp = [Sobol_temp;AmericanPutPrice_withoutIS,AmericanPutPrice_time,...
    AmericanPutPrice_paths];
%}
%%{
%% Generate LaTex code for a table of "call and put options",pathtype='QE',cubmethod='IID_MC'
clear input;
%temp = NaN(1,6);
call_QE = [IIDCallPrice_QE,IIDCallElapsed,1e6,1e6,1e6];%,SobolCallPrice_QE,temp];
put_QE = [IIDPutPrice_QE,IIDPutElapsed,1e6,1e6,1e6];%,SobolPutPrice_QE,temp];
content = [call_QE;iid_temp(1:2,:); put_QE;iid_temp(3:end,:)];
% content = [call_QE;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_QE;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'K=20','K=60','K=100','K=20','K=60','K=100','K=20','K=60','K=100'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'ECQE','ECw/IS','ECw/oIS','EPQE','EPw/IS',...
'EPw/oIS','APw/IS','APw/oIS'};

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
%}
%{
%% Generate LaTex code for a table of "call and put options",pathtype='QE',cubmethod='Sobol'
clear input;
% temp = NaN(1,6);
% call_QE = [IIDCallPrice_QE,temp];%,SobolCallPrice_QE,temp];
% put_QE = [IIDPutPrice_QE,temp];%,SobolPutPrice_QE,temp];
% content = [call_QE;Sobol_temp(1:2,:); put_QE;Sobol_temp(3:end,:)];
content = [Sobol_temp(1:2,:);Sobol_temp(3:end,:)];
% content = [call_QE;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_QE;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'K=20','K=60','K=100','K=20','K=60','K=100','K=20','K=60','K=100'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'ECw/IS','ECw/oIS','EPw/IS',...
'EPw/oIS','APw/IS','APw/oIS'};

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
%}
%% Reference

