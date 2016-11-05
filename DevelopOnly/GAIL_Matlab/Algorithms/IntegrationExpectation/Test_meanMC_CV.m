%% Initialization
% First we set up the basic common praramters for our examples.

gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for three months
inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0.05; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 120; %strike price
inp.payoffParam.putCallType = {'put'}; %looking at a put option
%inp.payoffParam.optType = {'amean'};
inp.priceParam.absTol = 0.02; %absolute tolerance of a two cents
inp.priceParam.relTol = 0; %zero relative tolerance
EuroPut = optPrice(inp); %construct an optPrice object 
disp('The price of the European put option')
disp(['    with a geometric Brownian motion is $' num2str(EuroPut.exactPrice,'%5.2f')])

%% The American put without control variates
% Next we create an American put |optPrice| object and use Monte Carlo to
% compute the price.

AmerPut = optPrice(EuroPut); %construct an American put object
AmerPut.payoffParam.optType = {'american'};
[AmerPutPrice,Aout] = genOptPrice(AmerPut)
disp(['The price of the American put option is $' ...
   num2str(AmerPutPrice,'%5.2f')])
disp(['   and this took ' num2str(Aout.time) ' seconds'])


%% The American put *with* control variates
% To use control variates we need to set up an |optPayoff| object with
% _two_ or more payoffs, the one whose expectation we want to compute, and the
% control variate(s)
AmerEuro = optPrice(inp);

AmerEuro.payoffParam = ...
   struct('optType',{{'american','euro'}}, ... %note two kinds of option payoffs
   'putCallType', {{'put','put'}}); %this needs to have the same dimension;
AmerEuro.priceParam.cubMethod = 'IID_MC_CV';
[AmerPutPrice,Aout] = genOptPrice(AmerEuro)



% [AmerEuroPrice, AEout] = meanMC_CV(@(n) (genOptPayoffs(AmerEuro,n)),...
%     AmerEuro.exactPrice(2:end),AmerEuro.priceParam.absTol,AmerEuro.priceParam.relTol,0.01); 
% disp(['The price of the American put option is $' ...
%    num2str(AmerEuroPrice,'%5.2f')])
% disp(['   and this took ' num2str(AEout.time) ' seconds,'])
% disp(['   which is ' num2str(AEout.time/Aout.time) ...
%    ' of the time without control variates'])
