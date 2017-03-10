%% This is a plot of Walsh coefficients in cone condition of reliable QMC
%% using arithematic mean Asian call option 

%% Initialization
% initialize the workspace and set option parameters
clc;clearvars;
inp.timeDim.timeVector = 1/52:1/52:4/52; %daily monitoring for 64 days 
inp.assetParam.initPrice = 120; %initial asset price
inp.assetParam.interest = 0.01; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 130; %strike price
inp.priceParam.absTol = 1e-3; %absolute tolerance
inp.priceParam.relTol = 0; %relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
EuroCall = optPrice(inp); %construct an optPrice object 
opt = optPayoff(EuroCall); %make a copy for no CV
% set the option type for target and CV functions 
    opt.payoffParam = struct( ...
    'optType',{{'amean'}},...
	'putCallType',{{'call'}}); 
% gather param for cubSobol_g 
abstol = inp.priceParam.absTol; % absolute tolerance
d = opt.timeDim.nSteps; 
% define payoff functions
f =@(x) genOptPayoffs(opt,x);

%% Checkings and plots
mmax = 14 
abstol = 1e-30;
m = 2^12;
l = 2^8;
l1 = l/2+1;
[~,~,y,kappanumap] =cubSobol_g(f,[zeros(1,d) ; ones(1,d)],'uniform',abstol,0, 'mmax',mmax);
figure
% set scale to log
set(gca,'xscale','log','yscale','log')
hold on
% plot all walsh coefficients(up to m), marker: black dots
scatter(1:m,abs(y(kappanumap(1:m))), '.','k');
% plot all points summed in 'diamond' , marker: green diamonds 
rangeDiamond = m:2^mmax;
s2 = scatter(rangeDiamond, abs(y(kappanumap(rangeDiamond))),'d', 'g');
% plot all points of 'square', marker: blue squres
rangeSquare = l1:l;
s3 = scatter(rangeSquare, abs(y(kappanumap(rangeSquare))), 's','b');
% plot all points summed in 'circle', marker: red circles 
rangeCircle = m:m:2^mmax;
s1 = scatter(rangeCircle, abs(y(kappanumap(rangeCircle))),'o','r');
hold off

% axis square
h_legend=legend([s1, s2, s3],...
       {'error bound',...
       '$\kappa \geq 2^{12}$',...
       '$\kappa \in [2^7, 2^8]$'},...
       'location','southwest');
set(h_legend,'FontSize',14);
set(h_legend,'interpreter','latex'); 
xlabel('$\kappa$','interpreter','latex')
ylabel('$|\hat{f}_{\kappa}|$','interpreter','latex')
legend('boxoff')
print('-deps', './cone.eps)
%
% %%
% _Author: Da Li 
