%% Estimation of the Expected Value Using |meanMC_CLT|
% For a random variable \(Y = f(\boldsymbol{X})\), we will estimate its expected value:
%
% \[
%  \mu=\mathbb{E}(Y) =
%  \mathbb{E}[f(\boldsymbol{X})]=\int_{\boldsymbol{R^d}} f(\boldsymbol{x}) \rho(\boldsymbol{x}) d\boldsymbol{x}
%  \approx \frac{1}{n}\sum_{i=1}^{n}{f(x_i)}, \, {x_i} \,
%  \text{IID} \sim
%  \, \rho
% \]
%
% We will approximate \(\mu\) using meanMC_CLT GAIL method. It is a
% IID Monte-Carlo algorithm using Central Limit Theorem. In order to improve
% computation efficiency, we will use control variates.
%%

function demoMCLT
%% Initialize the workspace and setting the display parameters
% These settings clean up the workspace and make the display beautiful.

gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters

%% Example 1: Estimate \(\mathbb{E}[f(\boldsymbol{X})]\) where \(f(\boldsymbol{x})=\exp(-\boldsymbol{x}^2)\) and \(\boldsymbol{X} \sim \mathcal{U} (0,1)\) using \(f(\boldsymbol{x})=\boldsymbol{x}\) as a control variate
absTol=1e-3; % absolute tolerance 
relTol=0; %relative tolerance
f=@(x)[exp(-x.^2),x];YXn=@(n)f(rand(n,1));%set up the random variable 

figure %plot f(x)
n=0:0.001:1;
plot(n,exp(-n.^2),'-'); 
ylabel('\(\exp(-x^2)\)')
xlabel('\(x\)')

s=struct('Y',YXn,'nY',1,'trueMuCV',1/2); % create a structure containing random variables, number of random variables and mean of the control variates 
[hmu,out]=meanMC_CLT(s,absTol,relTol);  % calculate the mean
exactsol=erf(1)*sqrt(pi)/2; %true mean
disp('Example 1')
disp(['Estimated mean is: ' num2str(hmu)])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nSample) ' points.'])
disp(['Real error was ' ...
    num2str(abs(exactsol-hmu))...
    ' which is less than the user input tolerance '...
    num2str(absTol) '.'])


%% Example 2: Price European Call option with stock price as a control variate
%Initialize option parameters for a European call option
inp.timeDim.timeVector = [1 2 3]; %%time increments
inp.assetParam.initPrice = 10; %initial stock price
inp.assetParam.interest = 0.01; %risk-free interest rate
inp.assetParam.volatility = 0.8; %volatility
inp.payoffParam.strike = 10; %strike price
inp.priceParam.absTol = 0.1; %absolute tolerance
inp.priceParam.relTol = 0; %relative tolerance
EuroCall = optPayoff(inp); % create a european call option
EuroCallPayoff=@(n) genOptPayoffs(EuroCall,n);
% Plot an empirical distribution of the European call option
n = 1e4; %number of payoffs to plot
payoffs = EuroCallPayoff(n); %generate n payoffs
sortedpay = sort(payoffs); %sort them
figure
plot(sortedpay,((1:n)-1/2)/n,'-'); %plot the empirical distribution function scenarios
xlabel('Payoff in dollars')
ylabel('CDF')
axis([0 300 0 1])
print -depsc PayoffCDF.eps %print the plot to a .eps file
% set up the stock price
load /Users/yueyili/GAIL_Dev/DevelopOnly/GAIL_Matlab/Algorithms/EnhancedAlgorithms/stockPriceHistory -ascii %load one year of stock price data into memory
S0 = stockPriceHistory(end); %stock price today
Delta = 1/250; %daily time increment in years
diffLogStockPrice = diff(log(stockPriceHistory)); %difference of the log of the stock prices
scDrift = mean(diffLogStockPrice); %sample mean
drift = scDrift/Delta; %estimated drift
scVolatility = std(diffLogStockPrice); %sample standard deviation
volatility = scVolatility/sqrt(Delta); %estimated volatility
timeFinal = 1/2; %final time
%estimate European call option
Y=@(n) [genOptPayoffs(EuroCall,n),S0*exp(drift*timeFinal+ volatility * sqrt(timeFinal)* randn(n,1))];
s=struct('Y',Y,'nY',1,'trueMuCV',S0);
[hmu,out]=meanMC_CLT(s,inp.priceParam.absTol,inp.priceParam.relTol);
disp('Example 2')
disp(['Estimated mean is: ' num2str(hmu)])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nSample) ' points.'])
disp(['Real error was ' ...
    num2str(abs(EuroCall.exactPrice-hmu))...
    ' which is less than the user input tolerance '...
    num2str(inp.priceParam.absTol) '.'])
%% Example 3: Keister's multidimensional integration
% We will evaluate the Keister's integration $I$ using meanMC_CLT.  Note
% that the we do a change of variable \(\boldsymbol{t} = \boldsymbol{x}/a\)
% and transform the integral:
%
% \begin{align*} I &= \int_{\mathbb{R}^d} \cos(a \lVert \boldsymbol{t}
% \rVert) \exp(-a^2 \lVert \boldsymbol{t} \rVert^2) \, a^d \mathrm{d}
% \boldsymbol{t}, \qquad a > 0, \\ & = \int_{\mathbb{R}^d}
% \underbrace{(2\pi a^2)^{d/2} \cos(a \lVert \boldsymbol{t} \rVert)
% \exp((1/2-a^2) \lVert \boldsymbol{t} \rVert^2)}_{f(\boldsymbol{t})}
% \times \underbrace{\frac{\exp(-\lVert \boldsymbol{t} \rVert^2/2)}
% {(2\pi)^{d/2}}}_{\varrho(\boldsymbol{t})} \, \mathrm{d} \boldsymbol{t} \\
% & = \mathbb{E}[f(\boldsymbol{T})], \qquad \text{where } \boldsymbol{T} \sim \mathcal{N}(\boldsymbol{0},
% \mathsf{I}). \end{align*}


%%
% define an anonymous function \(f\) as follows:
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
abstol = 0; %absolute tolerance
reltol = 0.01; %relative tolerance
dvec = 1:5; %vector of dimensions
avec = [1 1/sqrt(2)]; %default value of a 
IMCvec = zeros(size(dvec)); %vector of answers
f2= @(t,d) cell2mat(arrayfun(@(a) f(t,a,d),avec,'UniformOutput',false)); %a vector of funcion for each value of a
outT = zeros(size(dvec));%vector of time
outN=zeros(size(dvec));%vector of points
 for d = dvec
 f3=@(t)f2(t,d);%integration in dimension d
 YXn=@(n)f3(randn(n,d));%random generator
 s=struct('Y',YXn,'nY',size(avec,2)); 
 [IMCvec(d),out]= meanMC_CLT(s,abstol,reltol);
 outT(d)=out.time;
 outN(d)=out.nSample;
 end
[~,Ivec] = Keistertrue(dvec(end)); %true integration
relErrMC = abs(Ivec-IMCvec)./abs(Ivec);
 disp('Example 3')
disp(['The estimated integration for dimension ' num2str(dvec) ': ' num2str(IMCvec) ])
disp(['The algorithm took ' num2str(outT) ' seconds and '...
    num2str(outN) ' points.'])
disp(['Real error was ' ...
    num2str(relErrMC)...
    ' which is less than the user input tolerance '...
    num2str(reltol) '.'])


 
end
