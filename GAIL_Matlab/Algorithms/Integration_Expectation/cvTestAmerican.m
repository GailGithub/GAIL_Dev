%% Pricing Options with Control Variates
% This is a test for pricing american put option with CV. 

%% Initialization
% basic common praramters for our examples.
clc;clearvars;
iter=1; % number of trials
mmin=12; mmax=30;
fudge=@(m) 5*2.^-m;
update=1;% update beta or not

[da_out, da_out0, da_out1, da_out2, da_out3,...
 err, err0, err1, err2, err3]=deal([]); % container for outputs

% Initialize the workspace and the display parameters
inp.timeDim.timeVector = 1/52:1/52:16/52;%daily monitor for 64 days
inp.assetParam.initPrice = 36; %initial stock price
inp.assetParam.interest = 0.06; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 40; %strike price
inp.priceParam.absTol = 1e-3; %absolute tolerance of a nickel
inp.priceParam.relTol = 0; %zero relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
EuroCall = optPrice(inp); %construct an optPrice object 
opt = optPayoff(EuroCall); %make a copy
opt1 = optPayoff(EuroCall); %make a copy
opt2 = optPayoff(EuroCall); %make a copy

opt.payoffParam = struct( ...
	'optType',{{'american'}},...
	'putCallType',{{'put'}}); 
opt1.payoffParam = struct( ...
	'optType',{{'american','euro'}},...
	'putCallType',{{'put','put'}}); 
opt2.payoffParam = struct( ...
	'optType',{{'american','euro','stockprice'}},...
	'putCallType',{{'put','put',''}}); 

% make param shorter
abstol = inp.priceParam.absTol;
reltol = inp.priceParam.relTol;
d = opt.timeDim.nSteps; 
% get exact price of the option, 0 if not known

f =@(x) genOptPayoffs(opt,x);
f1.func =@(x) genOptPayoffs(opt1,x);
f1.cv = opt1.exactPrice;f1.cv=f1.cv(2:end); 
f2.func =@(x) genOptPayoffs(opt2,x);
f2.cv = opt2.exactPrice;f2.cv=f2.cv(2:end); 

%f3.func = {@(x) genOptPayoffs(opt,x),  @(x) genOptPayoffs(opt3,x)};
%f3.cv = [opt3.exactPrice];

fprintf(' abstol=%s \n ', abstol);
fprintf(' # of tests=%d \n ', iter);

% begin testing cubSobolcv no cv
for i=1:iter
	[q0,out0] = cubSobol_american_g(f,[zeros(1,d) ; ones(1,d)],abstol,0,mmin,mmax,'fudge',fudge);
    exactf = opt.exactPrice; exactf(isnan(exactf))=q0; 
    err0=[err0;abs(q0-exactf)];
    da_out0 = [da_out0;[out0.time, out0.n, out0.exitflag]];
end
fprintf('\n Results of cubSobolcv_g: \n');
fprintf('q0=%8.7f  \n',q0);
fprintf('max err(|q0-exact|): %.10f  \n', max(err0));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out0(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out0(:,2)));
fprintf('accumulated exitflag: %4d \n', sum(da_out0(:,3)));

% begin testing single cv
for i=1:iter
    [q1,out1] = cubSobol_american_g(f1,[zeros(1,d) ; ones(1,d)],abstol,0,mmin,mmax,'fudge',fudge,'betaUpdate',update);
    exactf = opt.exactPrice; exactf(isnan(exactf))=q0; 
    err1=[err1;abs(q1-exactf)];
    da_out1 = [da_out1;[out1.time, out1.n, out1.exitflag]];
end
fprintf('\n Results of cubSobolcv_g(single cv): \n');
fprintf('q1=%8.7f  \n',q1);
fprintf('max err(|q1-exact|): %.10f  \n', max(err1));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out1(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out1(:,2)));
fprintf('accumulated exitflag: %4d \n', sum(da_out1(:,3)));
fprintf('beta: %.4f \n', out1.beta);

% begin testing double cv
for i=1:iter
    [q2,out2] = cubSobol_american_g(f2,[zeros(1,d) ; ones(1,d)],abstol,0,mmin,mmax,'fudge',fudge,'betaUpdate',update);
    exactf = opt.exactPrice; exactf(isnan(exactf))=q0; 
    err2=[err2;abs(q2-exactf)];
    da_out2 = [da_out2;[out2.time, out2.n, out2.exitflag]];
end
fprintf('\n Results of cubSobolcv_g(double cv): \n');
fprintf('q2=%8.7f  \n',q2);
fprintf('max err(|q2-exact|): %.10f  \n', max(err2));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out2(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out2(:,2)));
fprintf('accumulated exitflag: %4d \n', sum(da_out2(:,3)));
fprintf('beta: %.4f \n', out2.beta);

%
% %%
% _Author: Da Li 
