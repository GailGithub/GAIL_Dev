%%{
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
inp.assetParam.pathType = 'QE_m'; 
inp.priceParam.cubMethod = 'IID_MC';
%Set optPayoff parameter
inp.assetParam.initPrice = 60;            % strike price
%Set error tolerance
inp.priceParam.absTol = 0;              % absolute tolerance
inp.priceParam.relTol = 0.01;           % three penny on the dollar relative tolerance
strike = [20,60,100];
avg = 5; % run five times and take average
%% European call option
%
% * cubMethod = 'IID_MC'

IIDCallPrice = zeros(1,3);
IIDCallPrice_paths = zeros(1,3);
IIDCallPrice_time = zeros(1,3);

for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourIIDCallPrice = optPrice(inp);
    for j = 1:avg        
        [temp, out] = genOptPrice(ourIIDCallPrice);
        IIDCallPrice(i) = IIDCallPrice(i) + temp;
        IIDCallPrice_paths(i) = IIDCallPrice_paths(i) + out.nPaths;
        IIDCallPrice_time(i) = IIDCallPrice_time(i) + out.time;
    end
end
IIDCallPrice = IIDCallPrice ./ avg;
IIDCallPrice_paths = IIDCallPrice_paths ./ avg;
IIDCallPrice_time = IIDCallPrice_time ./ avg;

iid_temp = [IIDCallPrice,IIDCallPrice_time,IIDCallPrice_paths];

inp.assetParam.pathType = 'QE'; 
IIDCallPrice= zeros(1,3);
IIDCallPrice_paths = zeros(1,3);
IIDCallPrice_time = zeros(1,3);

for i=1:3
    inp.payoffParam.strike = strike(i);
    %Construct an optPrice object
    ourIIDCallPrice = optPrice(inp);
    for j = 1:avg        
        [temp, out] = genOptPrice(ourIIDCallPrice);
        IIDCallPrice(i) = IIDCallPrice(i) + temp;
        IIDCallPrice_paths(i) = IIDCallPrice_paths(i) + out.nPaths;
        IIDCallPrice_time(i) = IIDCallPrice_time(i) + out.time;
    end
end
IIDCallPrice = IIDCallPrice ./ avg;
IIDCallPrice_paths = IIDCallPrice_paths ./ avg;
IIDCallPrice_time = IIDCallPrice_time ./ avg;

iid_temp = [IIDCallPrice,IIDCallPrice_time,IIDCallPrice_paths;iid_temp];


dif = NaN(1,9);
dif(1:3) = abs(iid_temp(1,1:3)-iid_temp(2,1:3));
%}
%% Generate LaTex code for a table of "call and put options",pathtype='QE',cubmethod='IID_MC'
clear input;
%temp = NaN(1,6);
content = [iid_temp;dif];
% content = [call_QE;iid_temp(1:2,:),Sobol_temp(1:2,:); ...
%     put_QE;iid_temp(3:4,:),Sobol_temp(3:4,:)];
input.data = content;
% Set column labels (use empty string for no label):
input.tableColLabels = {'K=20','K=60','K=100','K=20','K=60','K=100','K=20','K=60','K=100'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'Modified QE', 'Modified QE\_m','diff'};

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
input.tableCaption = 'Compare modified QE and modified QE with martingale correction';
% LaTex table label:
input.tableLabel = 'MyTableLabel';
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Switch to landscape table:
input.landscape = 0;
% call latexTable:
latex = latexTable(input);
