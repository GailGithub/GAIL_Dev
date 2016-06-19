%% Test Cases for class meanMC 

%% Example 1
% In this example, we are interested in estimating E[exp(U)], where U is
% uniformly distributed over [0, 1].
%
% The control variate used by CV is U.

num = 10; % produce 10 estimates
mu_cv = ones(num,1); 
mu = ones(num,1); 
time_cv = ones(num,1);
time = ones(num,1);


% setup for meanMC object
inp.in_param.abstol = 2e-4;
inp.in_param.reltol = 0;
inp.method = {'plain'};
inp.Yrand = @(n) exp(rand(n, 1));
inp.cv_param.YXrand = @expr;
inp.cv_param.muX = 0.5;

truemu = exp(1) - 1; % true answer

test1 = meanMC(inp);

for n = 1:num
    [tmu, param] = genMu(test1); 
    mu(n,1) = tmu;
    time(n,1) = param.time;  
end

test1.method = {'cv'}; % change method

for n = 1:num
    [tmu_cv, param_cv] = genMu(test1); 
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time; 
end


% plot estimates against running time for two methods

plot(mu_cv, time_cv, 'o', mu, time,'o') 
legend('cv','plain')
xlabel('estimates')
ylabel('time')
center = line([truemu truemu], [0 max(time)+0.5]);
center.Color = 'k';
left = line([truemu-test1.in_param.abstol  truemu-test1.in_param.abstol], ...
    [0 max([time' time_cv'])+0.5]);
right = line([truemu+test1.in_param.abstol  truemu+test1.in_param.abstol], ...
    [0 max([time' time_cv'])+0.5]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';
%% Example 2
% In this example, we are interested in estimating the average distance
% between two points in a 2-dimensional unit space using meanMCCV_g and
% meanMC_g.
%
% \[ \mu = \int_{[0,1]^2 \times [0,1]^2} \sqrt{(x_1-y_1)^2 +
% (x_2-y_2)^2} \, {\rm d}x_1 \, {\rm d}x_2 \, {\rm d}y_1 \, {\rm d}y_2\]
%
% The control variate used by meanMCCV_g is \( X_1, X_2, Y_1, Y_2 \).

abstol = 2e-4;

distfun1 = @distfun;
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

% plot estimates against running time for two methods

plot(mu_cv, time_cv, 'o', mu, time,'o') 
legend('meanMCCV\_g','meanMC\_g')
xlabel('estimates')
ylabel('time')
center = line([truemu truemu], [0 max(time_cv)+0.5]);
center.Color = 'k';
left = line([truemu-abstol  truemu-abstol], [0 max(time_cv)+0.5]);
right = line([truemu+abstol  truemu+abstol], [0 max(time_cv)+0.5]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

%% Example 3
% In this example, we are interested in pricing an Asian geometric call
% option with a strke price of 12. Its payoff depends on S(1), S(2) and
% S(3).
%
% The control variate used by meanMCCV_g is S(1), S(2) and S(3).

abstol = 0.01;

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

% plot estimates against running time for two methods

plot(mu_cv, time_cv, 'o', mu, time,'o') 
legend('meanMCCV\_g','meanMC\_g')
xlabel('estimates')
ylabel('time')
center = line([truemu truemu], [0 max(time)+0.5]);
center.Color = 'k';
left = line([truemu-abstol  truemu-abstol], [0 max(time)+0.5]);
right = line([truemu+abstol  truemu+abstol], [0 max(time)+0.5]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';


