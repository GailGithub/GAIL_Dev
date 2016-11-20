%% Initialization
% We start by setting up the basic common praramters for our examples.

gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for three months
inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0.05; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 120; %strike price
inp.payoffParam.putCallType = {'put'}; %looking at a put option
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
[AmerPutPrice,Aout] = genOptPrice(AmerPut);
disp(['The price of the American put option is $' ...
   num2str(AmerPutPrice,'%5.2f') ' +/- $' num2str(max(AmerPut.priceParam.absTol, ...
   AmerPut.priceParam.relTol*AmerPutPrice)) ])
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
[AmerEuroPutPrice,AEout] = genOptPrice(AmerEuro);

disp(['The price of the American put option is $' ...
   num2str(AmerEuroPutPrice,'%5.2f') ' +/- $' num2str(max(AmerEuro.priceParam.absTol, ...
   AmerEuro.priceParam.relTol*AmerEuroPutPrice)) ])
disp(['   and this took ' num2str(AEout.time) ' seconds,'])
disp(['   which is ' num2str(AEout.time/Aout.time * 100,'%5.2f') ...
   ' of the time without control variates'])


%% The Asian put without control variates
% Asian options, also known as arithmetic mean options, have payoff
% depending on the average of the stock price, not the final stock price.
% The typical discounted payoffs are as follows.
%
% \[
% \begin{array}{rcc}
% & \textbf{call} & \textbf{put} \\ \hline
% \textbf{payoff} & 
% \displaystyle \max\biggl(\frac 1d \sum_{j=1}^d S(jT/d) - K,0 \biggr)\mathsf{e}^{-rT} & 
% \displaystyle \max\biggl(K - \frac 1d \sum_{j=1}^d S(jT/d),0 \biggr)\mathsf{e}^{-rT} 
% \end{array}
% \]
AsianPut = optPrice(EuroPut);
AsianPut.payoffParam.optType = {'amean'};
[AsianPutPrice,AMeanout] = genOptPrice(AsianPut);
disp(['The price of this Asian arithmetic mean put option is $' num2str(AsianPutPrice,'%5.2f') ...
   ' +/- $' num2str(max(AsianPut.priceParam.absTol, AsianPut.priceParam.relTol*AsianPutPrice)) ])
disp(['   and it took ' num2str(AMeanout.time) ' seconds to compute']) %display results nicely


%% The Asian put *with* control variates
AGMean = optPrice(EuroPut);
AGMean.payoffParam = ...
    struct('optType',{{'amean','gmean'}},'putCallType',{{'put','put'}});
AGMean.priceParam.cubMethod = 'IID_MC_CV';
[AMeanPrice,AGout] = genOptPrice(AGMean);

disp(['The price of the Asian arithmetic mean put option is $' ...
   num2str(AMeanPrice,'%5.2f') ' +/- $' num2str(max(AGMean.priceParam.absTol, ...
   AGMean.priceParam.relTol*AMeanPrice)) ])
disp(['   and this took ' num2str(AGout.time) ' seconds,'])
disp(['   which is ' num2str(AGout.time/AMeanout.time * 100,'%5.2f') ...
   '% of the time without control variates'])


%% The Lookback call without control variates
% Lookback options consider the minimum or maximum asset price as their
% strike. The discounted payoffs are the following.
%
% \[
% \begin{array}{rcc}
% & \textbf{call} & \textbf{put} \\ \hline
% \textbf{payoff} & 
% \displaystyle \Bigl(S(T) - \min_{0 \le t \le T} S(t),0 \Bigr)\mathsf{e}^{-rT} & 
% \displaystyle \Bigl(\max_{0 \le t \le T} S(t) - S(T),0 \Bigr)\mathsf{e}^{-rT} 
% \end{array}
% \]

LookCall = optPrice(EuroPut);
LookCall.payoffParam = struct('putCallType',{{'call'}},'optType',{{'look'}});
[LookCallPrice,lout] = genOptPrice(LookCall);
disp(['The price of this Lookback call option is $' num2str(LookCallPrice,'%5.2f') ...
   ' +/- $' num2str(max(LookCall.priceParam.absTol, ...
   LookCall.priceParam.relTol*LookCallPrice)) ])
disp(['   and it took ' num2str(lout.time) ' seconds to compute']) %display results nicely


%% The Lookback call *with* control variates
LookCallCV = optPrice(EuroPut);
LookCallCV.payoffParam = ...
    struct('optType',{{'look','stockprice'}},'putCallType',{{'call','call'}});
LookCallCV.priceParam.cubMethod = 'IID_MC_CV';
[LookCallCVPrice,LCVout] = genOptPrice(LookCallCV);

disp(['The price of the Lookback call option is $' ...
   num2str(LookCallCVPrice,'%5.2f') ' +/- $' num2str(max(LookCallCV.priceParam.absTol, ...
   LookCallCV.priceParam.relTol*LookCallCVPrice)) ])
disp(['   and this took ' num2str(LCVout.time) ' seconds,'])
disp(['   which is ' num2str(LCVout.time/lout.time * 100,'%5.2f') ...
   '% of the time without control variates'])


%% Control variates for integral evaluation
% Suppose we want to evaluate $I = \int^2_0 \mathsf{e}^{-x^2}dx$, which is
% equivalent to $I = 2\int^1_0 \mathsf{e}^{-(2u)^2}du$. We then can create
% $U[0,1]$ IID samples to approximate this integral by simple Monte Carlo.

absTol = 0.001;
relTol = 0;
I = @(n) 2.*exp(-4.*rand(n,1).^2);
[I_MC,Iout] = meanMC_g(I,absTol,relTol);
disp(['The approximation of I is ' num2str(I_MC,'%6.3f') ' +/- ' ...
    num2str(max(absTol, I_MC*relTol)) ])
disp(['   and it took ' num2str(Iout.time) ' seconds to compute'])

%%
% Next we implement this by control variates method by taking the control
% variate as $X = \mathsf{e}^{-2u}$.
%%
% To ensure that we use exactly the same random numbers for $I$ and $X$, we
% enclose the data as a function below.
%%
% 
%   function YX = Test_CV_integral(n)
%   u = rand(n,1);
%   YX = [2.*exp(-4.*u.^2) exp(-2.*u)];
%   end
% 

[I_MC_CV, ICVout] = meanMC_CV(@(n) Test_CV_integral(n),(exp(2)-1)/(2*exp(2)),absTol,relTol);
disp(['The approximation of I is ' num2str(I_MC_CV,'%6.3f') ' +/- ' ...
    num2str(max(absTol, I_MC_CV*relTol)) ])
disp(['   and it took ' num2str(ICVout.time) ' seconds to compute'])
disp(['   which is ' num2str(ICVout.time/Iout.time * 100,'%5.2f') ...
   '% of the time without control variates'])
