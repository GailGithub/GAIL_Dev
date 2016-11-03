%% Pricing Digital Option with Monte Carlo and Quasi-Monte Carlo Methods
%
% An option, like a stock or bond, is a security. This type of investment 
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
% P, & S_{(T)} > K \\
% 0, & S_{(T)} \leq K\\
% \end{cases}\\
%  & \text{Put}\\
%  & \mathbb{1}_{[0,K]}(S_{(T)})=
% \begin{cases}
% P, & K > S_{(T)} \\
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
% To better explain, in the case of a Cash-or-nothing Call Option, if the 
% Stock Price at the maturity time \(S_{(T)}\) is above
% the Strike Price \(K\), the payoff of the option is _P_. 
% It does not depend if it closes $0.01 or $100.00 above the line, the 
% digital option is still worth the same amount of money. If \(S_{(T)}\) 
% closes below the underlying, then the option expires worthless, or 
% nothing at all.

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
% \Phi(x) & =\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{x}e^{-\frac{1}{2}z^2} \, 
% {\rm
% d}z \\
% \Phi(x) & = \text{Cumulative Distribution Function of the Normal 
% Distribution} \\ 
% d_1 &=\frac{\ln \Bigl(\frac{S_0}{K} \Bigr)+(r-q+\sigma^2/2)T}
% {\sigma\sqrt{T}} \\
% d_2 &=d_1-\sigma\sqrt{T} \\
% S_0 & = \text{initial stock price} \\
% K & = \text{strike price}\\
% T & = \text{maturity time}\\
% q & = \text{dividend rate}\\
% r & = \text{risk-free interest rate}\\
% \sigma & = \text{volatility}\\
% P & = \text{Payoff}\\
% \end{align*}
%
%% Pricing Digital Option for Cash-or-nothing
% Object Oriented Programming will be used in MATLAB. Every variable belongs
% to a class.
% _optPrice_ is a MATLAB(R) class that generates the prices for options. 
% _optPrice_ uses many different variables and it has many _superclasses_:
superclasses(optPrice)
%%
%
% In order to show how to price options, the first step is to create an 
% _input_ structure saying the parameters that we want to analyze.
%
% Creating the _input_ structure

clear; close all; clc;                     %Start Cleanning

%Payoff Parameters from optPayoff class
input.payoffParam.optType = {'digitalcash'};   %Cash-or-nothing Option Type
input.payoffParam.putCallType = {'call'};      %Call Option
input.payoffParam.strike = 12;                 %Strike Price ($)
input.payoffParam.digitalPay = 1;              %Payoff for one Option ($)

%Asset Path Parameters from assetPath class
input.assetParam.initPrice = 11;               %Initial Stock Price ($)
input.assetParam.interest = 0.01;              %Interest Rate (%)
input.assetParam.volatility = 0.5;             %Volatility (%)

%Option Price Parameters
input.priceParam.absTol = 0;                   %Absolute Tolerance
input.priceParam.relTol = 0.002;               %Relative Tolerance (%)

%Stochastic Process
input.timeDim.timeVector = 1;                  %Time to maturity (years)

%%
% _DigOption1_ is the _optPrice_ class

tic, DigOption1 = optPrice(input), toc
%%
% The function _genOptPrice_ is called from the IID Monte Carlo method for 
% Finance in <https://code.google.com/p/gail/ GAIL> Repository(meanMC_g).
% _mc_price_ is the price of the option, generated by IID Monte Carlo.

n = 5;  %Number of Samples
mc_price = zeros(1,n);
for i = 1:n
    tic,
    mc_price(i) = genOptPrice(DigOption1);
    toc
end
mc_price

disp(['Expect these answers to be within ' ...
   num2str(max(DigOption1.priceParam.absTol, ...
   DigOption1.priceParam.relTol * DigOption1.exactPrice)) ' of the true answer.'])

disp('The actual errors are')
fprintf(' %12.8f',abs(mc_price - DigOption1.exactPrice))

%%
% In order to compare the exact price from Monte Carlo and Quasi-Monte
% Carlo Methods, the 'Sobol' price parameter is used to construct 
% quasi-random points set.

SobolDigOption1 = optPrice(input);
SobolDigOption1.priceParam.cubMethod = 'Sobol';

sobol_price = zeros(1,n);
for i = 1:n
    tic,
    sobol_price(i) = genOptPrice(SobolDigOption1);
    toc
end
sobol_price

disp(['Expect these answers to be within ' ...
   num2str(max(SobolDigOption1.priceParam.absTol, ...
   SobolDigOption1.priceParam.relTol * SobolDigOption1.exactPrice)) ' of the true answer.'])

disp('The actual errors are')
fprintf(' %12.8f',abs(sobol_price - SobolDigOption1.exactPrice))
%%
%%% Plot the payoffs

plot(DigOption1,1e4);
%% Pricing Digital Option for Asset-or-nothing
%
% Copying the input structure from _DigOption1_. _DigOption2_ is the 
% _optPrice_ class

%%
DigOption2 = optPrice(input);                       %Copying all the parameters from DigOption1
DigOption2.payoffParam.optType = {'digitalasset'}   %Changing the Option Type to Asset-or-nothing
%%
% _mc_price2_ is the price of the option, generated by IID Monte Carlo.

mc_price2 = zeros(1,n);
for i = 1:n
    tic, 
    mc_price2(i) = genOptPrice(DigOption2);
    toc
end
mc_price2
disp(['Expect these answers to be within ' ...
   num2str(max(DigOption2.priceParam.absTol, ...
   DigOption2.priceParam.relTol * DigOption2.exactPrice)) ' of the true answer.'])

disp('The actual errors are')
fprintf(' %12.8f',abs(mc_price2 - DigOption2.exactPrice))
%%
% In order to compare the exact price from Monte Carlo and Quasi-Monte
% Carlo Methods, the 'Sobol' price parameter is used to construct 
% quasi-random points set.

SobolDigOption2 = optPrice(input);
SobolDigOption2.payoffParam.optType = {'digitalasset'};
SobolDigOption2.priceParam.cubMethod = 'Sobol';

sobol_price2 = zeros(1,n);
for i = 1:n
    tic, 
    sobol_price2(i) = genOptPrice(SobolDigOption2);
    toc
end
sobol_price2

disp(['Expect these answers to be within ' ...
   num2str(max(SobolDigOption2.priceParam.absTol, ...
   SobolDigOption2.priceParam.relTol * SobolDigOption2.exactPrice)) ' of the true answer.'])

disp('The actual errors are')
fprintf(' %12.8f',abs(sobol_price2 - SobolDigOption2.exactPrice))
%%
%%% Plot the payoffs

plot(DigOption2,1e4);
%% References
%
% "Binary Option." Wikipedia. Accessed June 5, 2015.
%
% "What Are Binary Options?" Binary Options. Accessed June 16, 2015.
%
% "Options Basics: What Are Options? | Investopedia." Investopedia. December 2, 2003. Accessed June 22, 2015.
%%
%
% Authors: Hartur Santi and Tianci Zhu