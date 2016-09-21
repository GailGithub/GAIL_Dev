%% This is a test of efficiency of cubSobol with bad~good CV 
% Problem: pricing up and in barrier call option
% Target function: up and in barrier call option payoff function
% CV function : European call option

%% Initialization
% initialize the workspace and set parameters
clc;clearvars;
iter=10; % number of tests
[da_out, da_out1, da_out2,...
 err, err1, err2]=deal([]); % initiate container for outputs

inp.timeDim.timeVector = 1/250:1/250:64/250; %daily monitoring for 64 days 
inp.assetParam.initPrice = 120; %initial asset price
inp.assetParam.interest = 0.01; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 130; %strike price
inp.payoffParam.barrier = 160; %barrier price, reset this to 150,140,130
inp.priceParam.absTol = 1e-3; %absolute tolerance
inp.priceParam.relTol = 0; %relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
EuroCall = optPrice(inp); %construct an optPrice object 
opt = optPayoff(EuroCall); %make a copy for no CV
opt1 = optPayoff(EuroCall); %make a copy for CV1

% set the option type for target and CV functions 
    opt.payoffParam = struct( ...
	'optType',{{'upin'}},...
	'putCallType',{{'call'}}); 
	opt1.payoffParam = struct( ...
	'optType',{{'upin','euro'}},...  
    'putCallType',{{'call','call'}}); 

% gather param for cubSobol_g 
abstol = inp.priceParam.absTol; % absolute tolerance
reltol = inp.priceParam.relTol;% relative tolerance
d = opt.timeDim.nSteps; 
fudge =@(m) 4*2.^-m;
mmin=10; mmax=30; % setup start  

% get exact price of the option
exactf = opt.exactPrice;

% define payoff functions
f =@(x) genOptPayoffs(opt,x);
f1.func =@(x) genOptPayoffs(opt1,x);
f1.cv = opt1.exactPrice;f1.cv=f1.cv(2:end); 

%% Tests start 
fprintf('abstol=%s \n ', abstol);
fprintf('# of tests=%d \n ', iter);

% begin a couple of test trials to get the environment stabilized 
for i=1:5
	cubSobol_g(f,[zeros(1,d) ; ones(1,d)],'uniform',abstol,0,mmin,mmax,fudge);
end

% begin testing no CV 
for i=1:iter
	[q,out] = cubSobol_g(f,[zeros(1,d) ; ones(1,d)],'uniform',abstol,0,mmin, mmax,fudge);
    err = [err;abs(q-exactf)];
    da_out = [da_out;[out.time, out.n]];
end
fprintf('\n Results of cubSobol_g: \n');
fprintf('q=%.10f  \n', q); % approximation from cubSobol_g
fprintf('exact=%.10f  \n', exactf); % exact results from formula 
fprintf('max err(|q-exact|): %.10f  \n', max(err)); % max error of all trials
fprintf('avg time of cubSobol_g: %.4f \n', mean(da_out(:,1))); % average time cost 
fprintf('avg n of cubSobol_g: %.2f \n', mean(da_out(:,2))); % average sample size

% begin testing single CV 
for i=1:iter
    [q1,out1] = cubSobol_g(f1,[zeros(1,d) ; ones(1,d)],'uniform',abstol,0,mmin,mmax,fudge);
    err1=[err1;abs(q1-exactf)];
    da_out1 = [da_out1;[out1.time, out1.n]];
end
fprintf('\n Results of cubSobol_g(single CV): \n');
fprintf('q1=%.10f  \n',q1);
fprintf('max err(|q1-exact|): %.10f  \n', max(err1));
fprintf('avg time of cubSobol_g: %.4f \n', mean(da_out1(:,1))); 
fprintf('avg n of cubSobol_g: %.2f \n', mean(da_out1(:,2)));
fprintf('beta: %.4f \n', out1.beta); % beta of last trial 

%
% %%
% _Author: Da Li 
