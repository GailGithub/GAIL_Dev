%% This is a test for multivariate normal probability
% Problem: compute multivariate normal probabilities
% Target function: density function of multivariate normal distribution  
% CV function #1: same function as target with diagonal covariance matrix 
% CV function #2: another variation of diagonal covariance matrix 

%% Initialization
% initialize the workspace and set parameters
clc;clearvars;
iter=10; % number of trials 
[da_out, da_out1, da_out2 ...
 err, err1, err2]=deal([]); % container for monitored outputs

s = 6;% dimension
abstol =1e-7;
fudge =@(m) 4*2.^-m;
mmin = 14; mmax=30;
sig = 0.1; mu = 0; % set sig=0.4, 0.7 for tests in thesis
% set up the hypercube
l = zeros(1,s); u = ones(1,s);
hyperbox = [l;u];
% set up the covariance matrix
cov = sig*ones(s,s);cov(1:s+1:s*s) = 1;

% define target function
f = @(x) exp(-0.5*sum(x.*(cov\x')',2))/(sqrt(prod(eig(cov))*(2*pi)^s));
% get approximated solution for target function 
[exactf, errf] = mvncdf(l,u,mu,cov);

% define CV functions
f1.func = @(x) [exp(-0.5*sum(x.*(cov\x')',2))/(sqrt(prod(eig(cov))*(2*pi)^s)), mvnpdf(x,mu,eye(s))];
f1.cv = mvncdf(l,u,mu,eye(s));
f2.func = @(x) [exp(-0.5*sum(x.*(cov\x')',2))/(sqrt(prod(eig(cov))*(2*pi)^s)), mvnpdf(x,mu,eye(s)), mvnpdf(x,mu,sig*eye(s))];
f2.cv = [mvncdf(l,u,mu,eye(s)),mvncdf(l,u,mu,sig*eye(s))];
% begin testing 
for i=1:iter
    [q,out] = cubSobol_g(f,hyperbox,'uniform',abstol,0,mmin,mmax,fudge);
    err = [err;abs(q-exactf)];
    da_out = [da_out;[out.time, out.n]];
end
fprintf('\n Results of cubSobolcv_g: \n');
fprintf('q=%.10f  \n', q);
%fprintf('exact=%.10f  \n', exactf);
fprintf('max err(|q-exact|): %.10f  \n', max(err));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out(:,2)));

% begin testing with single variate 
for i=1:iter
    [q1,out1] = cubSobol_g(f1,hyperbox,'uniform',abstol,0,mmin,mmax,fudge);
    err1 = [err1;abs(q1-exactf)];
    da_out1 = [da_out1;[out1.time, out1.n]];
end
fprintf('\n Results of cubSobolcv_g(single cv): \n');
fprintf('q1=%.10f  \n',q1);
fprintf('max err(|q1-exact|): %.10f  \n', max(err1));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out1(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out1(:,2)));
fprintf('beta: %.4f \n', out1.beta);

% begin testing with double variate 
for i=1:iter
    [q2,out2] = cubSobol_g(f2,hyperbox,'uniform',abstol,0,mmin,mmax,fudge);
    err2 = [err2;abs(q2-exactf)];
    da_out2 = [da_out2;[out2.time, out2.n]];
end
fprintf('\n Results of cubSobolcv_g(double cv): \n');
fprintf('q2=%.10f  \n',q2);
fprintf('max err(|q2-exact|): %.10f  \n', max(err2));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out2(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out2(:,2)));
fprintf('beta: %.4f \n', out2.beta);

%
% %%
% _Author: Da Li 
