%% Digital Option in Matlab
%
% An option, like a stock or bond, is a security. This type of investiment 
% is a contract that gives the buyer the right, but not the obligation, to
% buy or sell an asset at a specific price on a certain date.
%
% There are two different types of options, *calls* and *puts*:
%
% * A call gives the holder the right to buy an asset at a certain price 
% within a specific period of time. Buyers of calls hope that the stock 
% will increase substantially before the option expires.
% * A put gives the holder the right to sell an asset at a certain price 
% within a specific period of time. Buyers of puts hope that the price of 
% the stock will fall before the option expires.
%
%% Payoff of Digital Option
% In contrast to ordinary financial options that typically have a 
% continuous spectrum of payoff, _Digital Option_ is the type of option in 
% which the payoff can take only two possible results: the Payoff itself
% or nothing at all.
%
% There are two types of digital options: the *cash-or-nothing 
% option* and the *asset-or-nothing option*. The following characteristic
% function describes it in math terms:
%
% \[
% \begin{array}{rcc}
% \text{Cash-or-nothing} & \text{Call}\\
%  & \mathbb{1}_{[K,\infty)}(S_{(T)})=
% \begin{cases}
% Digital Pay, & S_{(T)} > K \\
% 0, & S_{(T)} \leq K\\
% \end{cases}\\
%  & \text{Put}\\
%  & \mathbb{1}_{[0,K]}(S_{(T)})=
% \begin{cases}
% Digital Pay, & K > S_{(T)} \\
% 0, & K \leq S_{(T)}\\
% \end{cases} \\
% \\
% \\
% \text{Asset-or-nothing} & \text{Call}\\
% & \mathbb{1}_{[K,\infty)}(S_{(T)})=
% \begin{cases}
% S_{(T)}, & S_{(T)} > K \\
% 0, & S_{(T)} \leq K\\
% \end{cases}\\
% & \text{Put}\\
% & \mathbb{1}_{[0,K]}(S_{(T)})=
% \begin{cases}
% S_{(T)}, & K > S_{(T)} \\
% 0, & K \leq S_{(T)}\\
% \end{cases} \\
% \end{array}
% \]
% 
% To better explain, in the case of a Cash-or-nothing Call Option, if the Stock Price at the maturity time \(S_{(T)}\) is above
% the Strike Price \(K\), the payoff of the option is _Digital Pay_. It does not
% depend if it closes $0.01 or $100.00 above the line, the digital option is 
% still worth the same amount of money. If \(S_{(T)}\) closes below the underlying, then the option
% expires worthless, or nothing at all.

%% Exact Price of Digital Option 
% When the underlying stock is assumed to follow the Black Scholes model or
% geometric Brownian motion model, there are formulas for the Exact Price
% of the digital option in terms of the standard Gaussian distribution
% function.  Let
%
% \[
% \begin{array}{rcc}
% & \text{Call}  & \text{Put} \\
% \text{Cash-or-nothing} & P.e^{-rT}\Phi(d_2) & P.e^{-rT}\Phi(-d_2) \\
% \text{Asset-or-nothing} & S_0.e^{-qT}\Phi(d_1) & S_0.e^{-qT}\Phi(-d_1)
% \end{array}
% \]
%
% under the following notation:
%
% \begin{align*}
% \Phi(x) & =\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{x}e^{-\frac{1}{2}z^2} \, {\rm
% d}z \\
% \Phi(x) & = \text{Cumulative Distribution Function of the Normal Distribution} \\ 
% d_1 &=\frac{\ln \Bigl(\frac{S_0}{K} \Bigr)+(r-q+\sigma^2/2)T}{\sigma\sqrt{T}} \\
% d_2 &=d_1-\sigma\sqrt{T} \\
% S_0 & = \text{initial stock price} \\
% K & = \text{strike price}\\
% T & = \text{maturity time}\\
% q & = \text{dividend rate}\\
% r & = \text{risk-free interest rate}\\
% \sigma & = \text{volatility}\\
% \end{align*}
%
%% Pricing Digital Option for Cash-or-nothing
% _optPrice_ is a MATLAB(R) class that generates the prices for options. 
% _optPrice_ uses many different variables and it has many superclasses:
superclasses(optPrice)
%%
%
% In order to show how to price an options, the first step is to create an 
% _input_ structure saying the parameters that we want to analyze.
%
% Creating the _input_ structure

%Payoff Parameters from optPayoff class
input.payoffParam.optType = {'digitalcash'};   %Cash-or-nothing Option Type
input.payoffParam.putCallType = {'call'};      %Call Option
input.payoffParam.strike = 12;                 %Strike Price is $12.00
input.payoffParam.digitalPay = 1;              %Each Digital Option will pay $1.00 if it ends in the money

%Asset Path Parameters from assetPath class
input.assetParam.initPrice = 11;               %Inicial Stock Price is $9.00
input.assetParam.interest = 0.01;              %Interest Rate is 1%
input.assetParam.volatility = 0.5;             %Volatility is 50%

%Option Price Parameters
input.priceParam.absTol = 0;                   %Do not use Absolute Tolerance
input.priceParam.relTol = 0.002;               %Relative Tolerance as 0.1%, or $0.02 in $10.00

%Stochastic Process
input.timeDim.timeVector = 1;                  %Time to maturity is one year and a single step

%%
% _DigOption1_ is the _optPrice_ class

tic, DigOption1 = optPrice(input), toc
%%
% From the above formulas, the Exact Price for _DigOption1_ is:
exactResult_DigOption1 = DigOption1.exactPrice
%%
% In order to compare the price from the formulas and the price
% from the Monte Carlo method for Finance in <https://code.google.com/p/gail/ GAIL> Repository(meanMC_g), the 
% function _genOptPrice_ is called.
% _mcResult_DigOption1_ is the price of the option, generated by IID Monte Carlo.

tic, mcResult_DigOption1 = genOptPrice(DigOption1), toc
%% 
% As it can be seen the _exactPrice_ for _DigOption1_ is the same value of
% the _IID_MC_ method with 0.2% as relative tolerance. 
%%
%%% Plot the Payoffs
plot(DigOption1,1e4);
%%
% In this plot can be seen that almost 65% of the times, the option will
% close worthless and the investor will lose all the money. But for almost
% 35% of the times, the investor will get $1.00 for each option.
% Is common to trade option with lots of 100 options, so the profit for the
% investor, if the option ends in the money would be:
%
% Investment = 100.$0.34 = $34.00
%
% Payoff = 100.$1.00 = $100.00
%
% Profit = $100.00 - $34.00 = $66.00
% 
%% Pricing Digital Option for Asset-or-nothing
%
% Copying the input structure from _DigOption1_

%Payoff Parameters from optPayoff class
%input.payoffParam.optType = {'digitalcash'};   %Cash-or-nothing Option Type
%input.payoffParam.putCallType = {'call'};      %Call Option
%input.payoffParam.strike = 12;                 %Strike Price is $12.00
%input.payoffParam.digitalPay = 1;              %Each Digital Option will
%pay $1.00 if it ends in the money

%Asset Path Parameters from assetPath class
%input.assetParam.initPrice = 11;               %Inicial Stock Price is $9.00
%input.assetParam.interest = 0.01;              %Interest Rate is 1%
%input.assetParam.volatility = 0.5;             %Volatility is 50%

%Option Price Parameters
%input.priceParam.absTol = 0;                   %Do not use Absolute Tolerance
%input.priceParam.relTol = 0.002;               %Relative Tolerance as 0.1%, or $0.02 in $10.00

%Stochastic Process
%input.timeDim.timeVector = 1;                  %Time to maturity is one year and a single step

%_DigOption2_ is the _optPrice_ class
%%
DigOption2 = DigOption1;                            %Copying all the parameters from DigOption1
DigOption2.payoffParam.optType = {'digitalasset'}   %Changing the Option Type to Asset-or-nothing
%%
% From the above formulas, the Exact Price for _DigOption2_ is:
exactResult_DigOption2 = DigOption2.exactPrice

%%
% _mcResult_DigOption2_ is the price of the option, generated by IID Monte Carlo.

tic, mcResult_DigOption2 = genOptPrice(DigOption2), toc
%%
%%% Plot the payoffs
plot(DigOption2,1e4);
%%
% 
%% References
%
% "Binary Option." Wikipedia. Accessed June 5, 2015.
%
% "What Are Binary Options?" Binary Options. Accessed June 16, 2015.
%
% "Options Basics: What Are Options? | Investopedia." Investopedia. December 2, 2003. Accessed June 22, 2015.




