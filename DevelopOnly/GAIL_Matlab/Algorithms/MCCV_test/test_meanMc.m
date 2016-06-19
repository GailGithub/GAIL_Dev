%% Test Cases for class meanMC 

%% Example 1
% In this example, we are interested in estimating E[exp(U)], where U is
% uniformly distributed over [0, 1].
%
% The control variate used by CV is U.

num = 10; % produce 10 estimates

mu = ones(num,1); 
mu_cv = ones(num,1);
mu_av = ones(num,1);
mu_all = ones(num,1);
time = ones(num,1);
time_cv = ones(num,1);
time_av = ones(num,1);
time_all = ones(num,1);

truemu = 1 - cos(1); % true answer

% setup for meanMC object
obj1.in_param.abstol = 2e-4;
obj1.in_param.reltol = 0;
obj1.method = {'plain'}; 
obj1.Yrand = @(n) sin(rand(n, 1));
test1 = meanMC(obj1);

for n = 1:num
    [tmu, param] = genMu(test1); 
    mu(n,1) = tmu;
    time(n,1) = param.time;  
end

obj1.method = {'cv'}; % control variate
obj1.cv_param.YXrand = @sinr_cv;
obj1.cv_param.muX = 0.5;
test1 = meanMC(obj1);

for n = 1:num
    [tmu_cv, param_cv] = genMu(test1); 
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time; 
end

obj1.method = {'av'}; % antithetic variate
obj1.av_param.YYrand = @sinr_av;
test1 = meanMC(obj1);

for n = 1:num
    [tmu_av, param_av] = genMu(test1); 
    mu_av(n,1) = tmu_av;
    time_av(n,1) = param_av.time; 
end

test1.method = {'plain','cv','av'}; % all methods

for n = 1:num
    [tmu_all, param_all] = genMu(test1); 
    mu_all(n,1) = tmu_all;
    time_all(n,1) = param_all.time; 
end

% plot estimates against running time for the three methods
plot(mu, time,'o', mu_cv, time_cv, 'o', mu_av, time_av, 'go')
legend('plain','cv','av')
xlabel('estimates')
ylabel('time')
tmax = max([time' time_cv' time_av'])*1.1;
center = line([truemu truemu], [0 tmax]);
center.Color = 'k';
left = line([truemu-test1.in_param.abstol  truemu-test1.in_param.abstol], ...
    [0 tmax]);
right = line([truemu+test1.in_param.abstol  truemu+test1.in_param.abstol], ...
    [0 tmax]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

% compare the best two methods with the combined one
plot(mu_cv, time_cv, 'o', mu_av, time_av, 'o', mu_all, time_all, 'go')
legend('cv','av','combined')
xlabel('estimates')
ylabel('time')
tmax = max([time_cv' time_av' time_all'])*1.1;
center = line([truemu truemu], [0 tmax]);
center.Color = 'k';
left = line([truemu-test1.in_param.abstol  truemu-test1.in_param.abstol], ...
    [0 tmax]);
right = line([truemu+test1.in_param.abstol  truemu+test1.in_param.abstol], ...
    [0 tmax]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';
%% Example 2
% In this example, we are interested in estimating the average distance
% between two points in a 2-dimensional unit space.
%
% \[ \mu = \int_{[0,1]^2 \times [0,1]^2} \sqrt{(x_1-y_1)^2 +
% (x_2-y_2)^2} \, {\rm d}x_1 \, {\rm d}x_2 \, {\rm d}y_1 \, {\rm d}y_2\]
%
% The control variate used by meanMC is \( X_1, X_2, Y_1, Y_2 \).

% use cubSobol_g to estimate the true mean 
distfunx = @(x) sqrt(sum((x(:,1:2)  - x(:,3:4)).^2,2));
truemu = cubSobol_g(distfunx,[zeros(1,4); ones(1,4)],'uniform',1e-5,0); 

% setup for meanMC object
obj2.in_param.abstol = 5e-4;
obj2.in_param.reltol = 0;
obj2.method = {'plain'};
obj2.Yrand = @(n) sqrt(sum((rand(n,2)  - rand(n,2)).^2,2));
test2 = meanMC(obj2);

for n = 1:num
    [tmu, param] = genMu(test2); 
    mu(n,1) = tmu;
    time(n,1) = param.time;  
end

obj2.method = {'cv'}; % control variate
obj2.cv_param.YXrand = @distfun_cv;
obj2.cv_param.muX = repmat(0.5,1,4);
test2 = meanMC(obj2);

for n = 1:num
    [tmu_cv, param_cv] = genMu(test2); 
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time; 
end

obj2.method = {'av'}; % antithetic variate
obj2.av_param.YYrand = @distfun_av;
test2 = meanMC(obj2);

for n = 1:num
    [tmu_av, param_av] = genMu(test2); 
    mu_av(n,1) = tmu_av;
    time_av(n,1) = param_av.time; 
end

test2.method = {'plain','cv','av'}; % all methods

for n = 1:num
    [tmu_all, param_all] = genMu(test2); 
    mu_all(n,1) = tmu_all;
    time_all(n,1) = param_all.time; 
end

% plot estimates against running time for the three methods
plot(mu, time,'o', mu_cv, time_cv, 'o', mu_av, time_av, 'go')
legend('plain','cv','av')
xlabel('estimates')
ylabel('time')
tmax = max([time' time_cv' time_av'])*1.1;
center = line([truemu truemu], [0 tmax]);
center.Color = 'k';
left = line([truemu-test2.in_param.abstol  truemu-test2.in_param.abstol], ...
    [0 tmax]);
right = line([truemu+test2.in_param.abstol  truemu+test2.in_param.abstol], ...
    [0 tmax]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

% compare the best two methods with the combined one
plot(mu, time, 'o', mu_all, time_all, 'o')
legend('plain','combined')
xlabel('estimates')
ylabel('time')
tmax = max([time' time_all'])*1.1;
center = line([truemu truemu], [0 tmax]);
center.Color = 'k';
left = line([truemu-test2.in_param.abstol  truemu-test2.in_param.abstol], ...
    [0 tmax]);
right = line([truemu+test2.in_param.abstol  truemu+test2.in_param.abstol], ...
    [0 tmax]);
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


