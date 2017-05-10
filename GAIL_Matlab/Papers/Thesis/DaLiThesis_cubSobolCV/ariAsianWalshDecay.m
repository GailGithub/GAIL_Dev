%% This is a plot of Walsh coefficients of arithematic mean Asian call option with and without CV (geomatric mean Asian call option) 

%% Initialization
% initialize the workspace and set option parameters
clc;clearvars;

inp.timeDim.timeVector = 1/250:1/250:64/250; %daily monitoring for 64 days 
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
opt1 = optPayoff(EuroCall); %make a copy for CV1
opt2 = optPayoff(EuroCall); %make a copy for CV2

% set the option type for target and CV functions 
    opt.payoffParam = struct( ...
	'optType',{{'amean'}},...
	'putCallType',{{'call'}}); 
	opt1.payoffParam = struct( ...
	'optType',{{'amean','gmean'}},...  
    'putCallType',{{'call','call'}}); 

% gather param for cubSobol_g 
abstol = inp.priceParam.absTol; % absolute tolerance
d = opt.timeDim.nSteps; 

% define payoff functions
f =@(x) genOptPayoffs(opt,x);
f1.func =@(x) genOptPayoffs(opt1,x);
f1.cv = opt1.exactPrice;f1.cv=f1.cv(2:end); 

%% Checkings and plots
m = 2^8;
mmax = 20;
abstol = 1e-30;
[~,~,y,kappanumap1] =cubSobol_g(f,[zeros(1,d) ; ones(1,d)],'uniform',abstol,0, 'mmax',mmax);
[~,~,y1,kappanumap2] =cubSobol_g(f1,[zeros(1,d) ; ones(1,d)],'uniform',abstol,0, 'mmax',mmax);
figure
scatter(1:size(kappanumap1,1)/m,abs(y(kappanumap1(1:end/m))),'red','MarkerFaceColor','red')
set(gca,'xscale','log','yscale','log')
hold on
scatter(1:size(kappanumap2,1)/m,abs(y1(kappanumap2(1:end/m))),'MarkerEdgeColor',[0 .5 0],'MarkerFaceColor',[0 .5 0])

hold off
% axis square
h_legend=legend('amean Asian Call Option',...
       'with gmean Asian Call Option as control variates',...
       'location','southwest')
%set(h_legend,'FontSize',14);
xlabel('$\kappa$','interpreter','latex')
ylabel('$|\hat{f}_{\kappa}|$','interpreter','latex')
legend('boxoff')

xmin = 4;
xmax = size(kappanumap1,1)/m*1.5;
ymin = min(min(abs(y(kappanumap1(1:end/m))), abs(y1(kappanumap2(1:end/m)))))/4;
ymax = max(max(abs(y(kappanumap1(1:end/m))), abs(y1(kappanumap2(1:end/m)))))*4;
axis([xmin xmax ymin ymax])
set(gca,'Xtick',10.^(0:floor(log10(xmax))))
set(gca,'Ytick',10.^(ceil(log10(ymin)):1:0))
print('-depsc', './ariAsianWalshDecay.eps')
%
% %%
% _Author: Da Li 
