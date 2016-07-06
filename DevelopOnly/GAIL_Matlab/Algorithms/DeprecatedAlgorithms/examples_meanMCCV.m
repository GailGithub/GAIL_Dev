%% Examples using meanMCCV_g  
% Author: Tianpei Qian
%% Example 1
% In this example, we are interested in estimating E[exp(U)], where U is
% uniformly distributed over [0, 1].
%
% The control variate used by meanMCCV_g is U.

num = 10; % produce 10 estimates
mu_cv = ones(num,1); % a vector containing the estimated means using meanMCCV_g
mu = ones(num,1);  % a vector containing the estimated means using meanMC_g
time_cv = ones(num,1); % a vector containing time used by meanMCCV_g
time = ones(num,1); % a vector containing time used by meanMC_g

abstol = 2e-4; % absolute tolerance rate

expr1 = @expr_cv; % expr_cv is defined in the folder test_fns
expr2 = @(n) exp(rand(n, 1));

truemu = exp(1) - 1; % true answer

for n = 1:num
    [tmu_cv, param_cv] = meanMCCV_g(expr1, 0.5, expr2, abstol, 0); % use meanMCCV_g
    [tmu, param] = meanMC_g(expr2, abstol, 0); % use meanMC_g
    
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time;
    mu(n,1) = tmu;
    time(n,1) = param.time;   
end


% Let's plot estimates against running time for two methods 

plot(mu_cv, time_cv, 'o', mu, time,'o') 
legend('meanMCCV\_g','meanMC\_g')
xlabel('estimates')
ylabel('time')
center = line([truemu truemu], [0 max(time)*1.1]);
center.Color = 'k';
left = line([truemu-abstol  truemu-abstol], [0 max(time)*1.1]);
right = line([truemu+abstol  truemu+abstol], [0 max(time)*1.1]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

%%
% It can be seen that meanMCCV_g is much more efficient than meanMC_g.
%% Example 2
% In this example, we are interested in estimating the average distance
% between two points in a 2-dimensional unit space using meanMCCV_g and
% meanMC_g.
%
% \[ \mu = \int_{[0,1]^2 \times [0,1]^2} \sqrt{(x_1-y_1)^2 +
% (x_2-y_2)^2} \, {\rm d}x_1 \, {\rm d}x_2 \, {\rm d}y_1 \, {\rm d}y_2\]
%
% The control variate used by meanMCCV_g is \( X_1, X_2, Y_1, Y_2 \).

abstol = 2e-4; % absolute tolerance

distfun1 = @distfun_cv; % distfun_cv is defined in the folder test_fns
distfun2 = @(n) sqrt(sum((rand(n,2)  - rand(n,2)).^2,2));

distfunx = @(x) sqrt(sum((x(:,1:2)  - x(:,3:4)).^2,2));
truemu = cubSobol_g(distfunx,[zeros(1,4); ones(1,4)],'uniform',1e-5,0); 
% use cubSobol_g to estimate the true mean 

for n = 1:num   
    [tmu_cv, param_cv] = meanMCCV_g(distfun1,repmat(0.5,1,4), distfun2, abstol, 0); % use meanMCCV_g
    [tmu, param] = meanMC_g(distfun2,abstol,0); % use meanMC_g
    
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time;
    mu(n,1) = tmu;
    time(n,1) = param.time;   
end

% Let's plot estimates against running time for two methods
plot(mu_cv, time_cv, 'o', mu, time,'o') 
legend('meanMCCV\_g','meanMC\_g')
xlabel('estimates')
ylabel('time')
center = line([truemu truemu], [0 max(time_cv)*1.1]);
center.Color = 'k';
left = line([truemu-abstol  truemu-abstol], [0 max(time_cv)*1.1]);
right = line([truemu+abstol  truemu+abstol], [0 max(time_cv)*1.1]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

%%
% It can be seen that the two methods have almost the same performance. 
% This is because the control variates are almost uncorrelated with the response.
%% Example 3
% In this example, we are interested in pricing an Asian geometric call
% option with a strke price of 12. Its payoff depends on S(1), S(2) and
% S(3).
%
% The control variate used by meanMCCV_g is S(1), S(2) and S(3).

abstol = 0.01; % absolute tolerance

inp.payoffParam.optType = {'gmean'};
inp.payoffParam.putCallType = {'call'};
inp.payoffParam.strike = 12;
inp.assetParam.initPrice = 11;
GMean = optPrice(inp);
truemu = GMean.exactPrice; % true answer

r = 0.01; % interest rate
initP = 11; % initial price of asset
muX = [exp(r*1)*initP exp(r*2)*initP exp(r*3)*initP];

for n = 1:num   
    [tmu_cv, param_cv] = meanMCCV_g(@gmean_cv,muX,@gmean,abstol, 0); % use meanMCCV_g
    [tmu, param] = meanMC_g(@gmean,abstol,0); % use meanMC_g
    
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time;
    mu(n,1) = tmu;
    time(n,1) = param.time;   
end

% Let's plot estimates against running time for two methods

plot(mu_cv, time_cv, 'o', mu, time,'o') 
legend('meanMCCV\_g','meanMC\_g')
xlabel('estimates')
ylabel('time')
center = line([truemu truemu], [0 max(time)*1.1]);
center.Color = 'k';
left = line([truemu-abstol  truemu-abstol], [0 max(time)*1.1]);
right = line([truemu+abstol  truemu+abstol], [0 max(time)*1.1]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

%%
% It can be seen that meanMCCV_g is much more efficient than meanMC_g. 
