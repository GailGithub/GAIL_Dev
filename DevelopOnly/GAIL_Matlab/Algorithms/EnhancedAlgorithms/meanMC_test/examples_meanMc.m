%% Examples for class meanMC 
% Author: Tianpei Qian
%% Example 1
% In this example, we are interested in estimating E[sin(U)], where U is
% uniformly distributed over [0, 1].
%
% The control variate used by CV is U.

num = 10; % produce 10 estimates

mu = ones(num,1); % a vector containing the estimated means with method 'plain'
mu_cv = ones(num,1); % a vector containing the estimated means with method 'cv'
mu_av = ones(num,1); % a vector containing the estimated means with method 'av'
mu_all = ones(num,1); % a vector containing the estimated means with method 'combined'
time = ones(num,1); % a vector containing time used method 'plain'
time_cv = ones(num,1); % a vector containing time used method 'cv'
time_av = ones(num,1); % a vector containing time used method 'av'
time_all = ones(num,1); % a vector containing time used method 'combined'

sinr_help = @(x) [sin(x(:,1)) x(:,1)];
sinr_cv = @(n) sinr_help(repmat(rand(n,1),1,2));
sinr_av = @(n) sin(bsxfun(@plus, bsxfun(@times, rand(n, 1), [1, -1]), [0, 1]));

truemu = 1 - cos(1); % true answer

% setup for meanMC object
obj1.in_param.abstol = 2e-4; % absolute tolerance rate
obj1.in_param.reltol = 0;
obj1.method = {'plain'}; 
obj1.Yrand = @(n) sin(rand(n, 1));
test1 = meanMC(obj1);

genMu(test1); % run once 
for n = 1:num
    [tmu, param] = genMu(test1); 
    mu(n,1) = tmu;
    time(n,1) = param.time;  
end

obj1.method = {'cv'}; % control variate
obj1.cv_param.YXrand = sinr_cv; 
obj1.cv_param.muX = 0.5;
test1 = meanMC(obj1);

genMu(test1); % run once 
for n = 1:num
    [tmu_cv, param_cv] = genMu(test1); 
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time; 
end

obj1.method = {'av'}; % antithetic variate
obj1.av_param.YYrand = sinr_av; 
test1 = meanMC(obj1);

genMu(test1); % run once 
for n = 1:num
    [tmu_av, param_av] = genMu(test1); 
    mu_av(n,1) = tmu_av;
    time_av(n,1) = param_av.time; 
end

test1.method = {'plain','cv','av'}; % all methods

genMu(test1); % run once 
for n = 1:num
    [tmu_all, param_all] = genMu(test1); 
    mu_all(n,1) = tmu_all;
    time_all(n,1) = param_all.time; 
end

% plot estimates against running time for the three methods
figure()
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
figure()
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
%%
% In this example, the method 'cv' is the most efficient method. 
% The method 'combined' is only a little worse than the optimal method.
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
obj2.in_param.abstol = 2e-4;
obj2.in_param.reltol = 0;
obj2.method = {'plain'};
obj2.Yrand = @(n) sqrt(sum((rand(n,2)  - rand(n,2)).^2,2));
test2 = meanMC(obj2);

genMu(test2); % run once 
for n = 1:num
    [tmu, param] = genMu(test2); 
    mu(n,1) = tmu;
    time(n,1) = param.time;  
end

obj2.method = {'cv'}; % control variate
obj2.cv_param.YXrand = @distfun_cv; % dist_cv is defined in the folder test_fns
obj2.cv_param.muX = repmat(0.5,1,4);
test2 = meanMC(obj2);

genMu(test2); % run once 
for n = 1:num
    [tmu_cv, param_cv] = genMu(test2); 
    mu_cv(n,1) = tmu_cv;
    time_cv(n,1) = param_cv.time; 
end

obj2.method = {'av'}; % antithetic variate
obj2.av_param.YYrand = @distfun_av; % dist_av is defined in the folder test_fns
test2 = meanMC(obj2);

genMu(test2); % run once 
for n = 1:num
    [tmu_av, param_av] = genMu(test2); 
    mu_av(n,1) = tmu_av;
    time_av(n,1) = param_av.time; 
end

test2.method = {'plain','cv','av'}; % all methods

genMu(test2); % run once 
for n = 1:num
    [tmu_all, param_all] = genMu(test2); 
    mu_all(n,1) = tmu_all;
    time_all(n,1) = param_all.time; 
end

% plot estimates against running time for the three methods
figure()
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

% compare the best method with the combined one
figure()
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

%%
% In this example, the method 'plain' is the most efficient method. 
% Again, the method 'combined' is only a little worse than the optimal method.
%% Example 3
% In this example, we are interested in pricing an up-and-in call option
% option with a strke price of 12. Its payoff depends on S(1), S(2) and
% S(3).
% The control variates used are the payoffs of a European call options with the
% different strike prices.

% setup for meanMC object
obj3.in_param.abstol = 1e-1; % absolute error tolerance
obj3.in_param.reltol = 0; % relative error tolerance
obj3.method = {'plain'}; 

inp.payoffParam.optType = {'upin'};
inp.payoffParam.putCallType = {'call'};
inp.payoffParam.strike = 110;
inp.payoffParam.barrier = 180;
inp.assetParam.initPrice = 100;

upin = optPayoff(inp);

% use cubSobol to estimate truemu
inp.priceParam.cubMethod = {'Sobol'};
inp.priceParam.absTol = 5e-3;
upin_sobol = optPrice(inp);
truemu = genOptPrice(upin_sobol);

obj3.Yrand = @(n) genOptPayoffs(upin, n);
test3 = meanMC(obj3);

genMu(test3); % run once 
for n = 1:num
    [tmu, param] = genMu(test3); 
    mu(n,1) = tmu;
    time(n,1) = param.time;  
end

mu_cv1 = ones(num,1);
mu_cv2 = ones(num,1);
mu_cv3 = ones(num,1);
time_cv1 = ones(num,1);
time_cv2 = ones(num,1);
time_cv3 = ones(num,1);

% single control variate
obj3.method = {'cv'}; 
inp.payoffParam.optType = {'upin', 'euro'};
inp.payoffParam.putCallType = {'call', 'call'};

upin_cv1 = optPayoff(inp);

obj3.cv_param.YXrand = @(n) genOptPayoffs(upin_cv1, n);
obj3.cv_param.muX = upin_cv1.exactPrice(2);
test3 = meanMC(obj3);

genMu(test3); % run once 
for n = 1:num
    [tmu_cv1, param_cv1] = genMu(test3); 
    mu_cv1(n,1) = tmu_cv1;
    time_cv1(n,1) = param_cv1.time; 
end

% add multiple control variates
strikes = 110:10:180; 
nstrikes = length(strikes);
inp.payoffParam.optType = [{'upin'} repmat({'euro'},1,nstrikes)];
inp.payoffParam.putCallType = repmat({'call'},1,nstrikes+1);
inp.payoffParam.strike = [110 strikes];

upin_cv2 = optPayoff(inp);

obj3.cv_param.YXrand = @(n) genOptPayoffs(upin_cv2, n);
obj3.cv_param.muX = upin_cv2.exactPrice(2:(nstrikes+1));
test3 = meanMC(obj3);

genMu(test3); % run once 
for n = 1:num
    [tmu_cv2, param_cv2] = genMu(test3); 
    mu_cv2(n,1) = tmu_cv2;
    time_cv2(n,1) = param_cv2.time; 
end

% use Ridge regression to estimate control variate coefficients
obj3.cv_param.ridge = 0.2;
test3 = meanMC(obj3);

genMu(test3); % run once 
for n = 1:num
    [tmu_cv3, param_cv3] = genMu(test3); 
    mu_cv3(n,1) = tmu_cv3;
    time_cv3(n,1) = param_cv3.time; 
end

% plot estimates against running time for the plain method and single cv
figure()
plot(mu, time,'o', mu_cv1, time_cv1, 'o')
legend('plain','cv (1 control)')
xlabel('estimates')
ylabel('time')
tmax = max([time' time_cv1'])*1.1;
center = line([truemu truemu], [0 tmax]);
center.Color = 'k';
left = line([truemu-test3.in_param.abstol  truemu-test3.in_param.abstol], ...
    [0 tmax]);
right = line([truemu+test3.in_param.abstol  truemu+test3.in_param.abstol], ...
    [0 tmax]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';
%%
% In this case, the method 'cv' is much more efficient than the method
% 'plain'.

% compare different versions of cv
figure()
plot(mu_cv1, time_cv1, 'o', mu_cv2, time_cv2, 'o', mu_cv3, time_cv3, 'go')
legend('cv (1 control)','cv (8 controls)','cvRidge (8 controls)')
xlabel('estimates')
ylabel('time')
tmax = max([time_cv1' time_cv2' time_cv3'])*1.1;
center = line([truemu truemu], [0 tmax]);
center.Color = 'k';
left = line([truemu-test3.in_param.abstol  truemu-test3.in_param.abstol], ...
    [0 tmax]);
right = line([truemu+test3.in_param.abstol  truemu+test3.in_param.abstol], ...
    [0 tmax]);
left.Color = 'k';
left.LineStyle = ':';
right.Color = 'k';
right.LineStyle = ':';

%%
% It can be seen that adding multiple controls makes the 'cv' method more
% efficient. Also, using Ridge regression to estimate control vocariates
% does not improve the efficiency. This may be because we
% already have enough samples to estimate the coefficients.
