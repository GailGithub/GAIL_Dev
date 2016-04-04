%% demo_normal_probabilities_CSC
function demo_normal_probabilities
%% Computing Normal probabilities with covariance matrix cov
% Dimension of the problem. In the literature it sometimes appears as s.
d = 10; 

%% Integration parameters
% Set maximum budget permitted in the algorithm in the form 2^mmax. It
% depends on your computer memory. Also, set absolute error tolerance and
% relative error tolerance.
mmax = 22; 
abstol = 1e-3; 
reltol = 0;  % If 0, we have pure absolute error.

%% Dumb test with cov = Id 
% To check real error:
mu = zeros(d,1);
cov = eye(d);
factor = 3.5;
% Integration limits
hyperbox = [-factor*ones(1,d) ; factor*ones(1,d)]; 
% Solution approx_prob and integration output parameters in out_param
[approx_prob1,out_param1] = multi_normcdf(hyperbox,mu,cov,abstol,reltol,mmax); 
disp('Error checking cubSobol_g example 1')
disp(['Solution for d=' num2str(d) ', cubSobol_g : ' num2str(approx_prob1)])
disp(['Real error was ' ...
    num2str(abs((gail.stdnormcdf(factor)-gail.stdnormcdf(-factor))^d-approx_prob1))...
    ' which is less than the tolerance ' num2str(abstol) '. Time was '...
    num2str(out_param1.time) ' seconds.'])

%% Another calibration example 
% To check real error:
mu = zeros(d,1);
% off-diagional value of the covariance matrix
sig = 0.6; 
% definition of the covariance matrix
cov = sig*ones(d,d); 
% definition of the covariance matrix
cov(1:d+1:d*d) = 1; 
% Integration limits
hyperbox = [-Inf*ones(1,d) ; sqrt(d)*rand(1,d)];
% Solution approx_prob and integration output parameters in out_param
[approx_prob2,out_param2] = multi_normcdf(hyperbox,mu,cov,abstol,reltol,mmax);
% Computing the real solution
[exact_p2 , e_out2] = cubSobol_g(...
  @(t) prod(gail.stdnormcdf(bsxfun(@plus,hyperbox(2,:),sqrt(sig)*t)/sqrt(1-sig)),2),...
  [-Inf;Inf],'normal',abstol/10^3,0); 
disp('Error checking cubSobol_g example 2')
disp(['Solution for d=' num2str(d) ', cubSobol_g : ' num2str(approx_prob2)])
disp(['Real error was ' num2str(abs(exact_p2-approx_prob2)) ...
    ' which is less than the tolerance ' num2str(abstol) '. Time was '...
    num2str(out_param2.time) ' seconds.'])

%% Real computation
mu = zeros(d,1);
sig = 0.6; % off-diagional value of the covariance matrix
cov = sig*ones(d,d); % definition of the covariance matrix
cov(1:d+1:d*d) = 1; % definition of the covariance matrix

hyperbox = [-(d/3)*rand(1,d) ; (d/3)*rand(1,d)]; % Integration limits
% Solution approx_prob and integration output parameters in out_param
[approx_prob3,out_param3] = multi_normcdf(hyperbox,mu,cov,abstol,reltol,mmax); 
disp('Real example cubSobol_g')
disp(['Solution for d=' num2str(d) ', cubSobol_g : ' num2str(approx_prob3)])
disp(['Time was ' num2str(out_param3.time) ' seconds.'])

%% Testing with Monte Carlo
C = chol(cov)';
a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1); s = gail.stdnormcdf(a); 
e = gail.stdnormcdf(b);
Y = @(n) f(s,e,hyperbox,rand(n,d-1),C);
[approx_prob4,out_param4] = meanMC_g(Y,abstol,reltol,'tbudget',5000);
disp('Real example meanMC_g')
disp(['Solution for d=' num2str(d) ', meanMC_g : ' num2str(approx_prob4)])
disp(['Time was ' num2str(out_param4.time) ' seconds.'])
end


function [p,out, y, kappanumap] = multi_normcdf(hyperbox,mu,cov,abstol,reltol,mmax)
% multi_normcdf computes the cumulative distribution function of the
% multivariate normal with mean mu, covariance matrix cov and within the
% region defined by hyperbox.

hyperbox = bsxfun(@minus, hyperbox,mu');
C = chol(cov)'; d = size(C,1);
a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1); s = gail.stdnormcdf(a); 
e = gail.stdnormcdf(b);
[p, out, y, kappanumap] = cubSobol_g(...
    @(x) f(s,e,hyperbox,x,C), [zeros(1,d-1);ones(1,d-1)],...
    'uniform',abstol,reltol,'mmax',mmax);

end

function f_eval = f(s,e,hyperbox,w,C)
    f_eval = (e-s)*ones(size(w,1),1);
    aux = ones(size(w,1),1);
    y = [];
    for i = 2:size(hyperbox,2);
        y = [y gail.stdnorminv(s+w(:,i-1).*(e-s))];
        aux = sum(bsxfun(@times,C(i,1:i-1),y),2);
        a = (hyperbox(1,i)-aux)/C(i,i);
        b = (hyperbox(2,i)-aux)/C(i,i);
        s = gail.stdnormcdf(a);
        e = gail.stdnormcdf(b);
        f_eval = f_eval .* (e-s);
    end
end
