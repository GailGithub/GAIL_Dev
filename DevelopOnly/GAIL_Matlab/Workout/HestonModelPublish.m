%% Heston Stochastic Volatility Modeling
% 
% Stochastic volatility models play an important role in improving the 
% accuracy of pricing financial derivatives. There are several well-known 
% stochastic volatility models: the Hull-White model (1987), the 
% Scott-Chesny model (1989), the Heston model (1993) and the SABR model 
% (2002). The Heston model is of particular interest to us since it is one 
% of the most widely used stochastic volatility models.
%
% The Heston Model is defined by the following two stochastic differential
% equations:
%
% $d X(t)=\sqrt{V(t)}X(t)d W_X(t)$
%
% $dV(t)=\kappa (\theta-V(t))dt+\nu\sqrt{V(t)}dW_V(t)$
%
% where
%
% $X(t)=$ asset price process,
%
% $V(t)=$ instantaneous variance of relative changes to $X(t)$,
%
% $d W_1$ and $d W_2$ are two standard Brownian motions which are correlated, $\rho$ being the correlation,
%
% $\kappa$ is the speed of mean reversion,
%    
% $\theta$ is the value of the long-term variance,
%
% $\nu$  is the volatility of volatility.

%% Advantage of the Algorithm
% Our modified algorithm can be used for any $\nu$, the volatility of the
% asset price volatility, greater than or equal to zero. Besides, the
% modified algorithm can reduce the numerical error even when $\nu\approx0$

%% Generate Option prices using Quadratice Exponential Scheme
%Set assetPath parameters
T=1;                                    % end time
delta_t=0.1;                            % time increment
t0 = delta_t;                           % start time
inp.timeDim.timeVector = t0:delta_t:T;  % time vector
inp.assetParam.initPrice = 100;         % initial asset price
inp.assetParam.interest = 0.04;            % risk-free interest rate
inp.assetParam.volatility = 0.3;        % fixed vlatility of asset prices
inp.assetParam.Vinst = 0.09;            % initial value of volatility
inp.assetParam.Vlong = 0.09;            % theta
inp.assetParam.kappa = 1;               % kappa
inp.assetParam.nu = 0;                  % volatility of asset price volatility
inp.assetParam.rho = 0.5;               % rho

%Set optPayoff parameter
inp.payoffParam.strike = 90;            % strike price

%Set error tolerance
inp.priceParam.absTol = 0;              %absolute tolerance
inp.priceParam.relTol = 0.01;           %one penny on the dollar relative tolerance

%%
inp.assetParam.pathType = 'QE';         % path type QE
%Construct an optPrice object
ourQEPrice = optPrice(inp) 
QEPrice=zeros(1,5);
%Generate Option Price
for i=1:5
    tic,
    QEPrice(i) = genOptPrice(ourQEPrice);
    toc
end
QEPrice
% Calculate relative tolerance when nu=0
variances = 0;
if inp.assetParam.nu==0
    reldiff = abs(QEPrice-ourQEPrice.exactPrice)/ourQEPrice.exactPrice
end
variaces = sum((QEPrice-ourQEPrice.exactPrice).^2)/5

%% Generate Option prices using Quadratice Exponential Scheme with Martingale Correction
inp.assetParam.assetPath='QE_m';
ourQEmPrice = optPrice(inp) 
QEmPrice=zeros(1,5);
%Generate Option Price
for i=1:5
    tic,
    QEmPrice(i) = genOptPrice(ourQEmPrice);
    toc
end
QEmPrice
% Calculate relative tolerance when nu=0
if inp.assetParam.nu==0
    reldiffm = abs(QEmPrice-ourQEmPrice.exactPrice)/ourQEmPrice.exactPrice
end
variaces = sum((QEmPrice-ourQEmPrice.exactPrice).^2)/5
%% Reference
%
% Andersen, Leif B. G. "Efficient Simulation of the Heston Stochastic Volatility Model."
