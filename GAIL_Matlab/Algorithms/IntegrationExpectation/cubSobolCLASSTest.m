% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:
% 
w.f = @(x) prod(x,2); w.domain = [zeros(1,2);ones(1,2)]; 
w.absTol=1e-5; w.relTol=0;
q = cubSobol_gCLASS(w);

% Example 2: NOT OK
% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2
% in the interval R^3 where x1, x2 and x3 are normally distributed:
% 
w.f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; w.domain = [-inf(1,3);inf(1,3)];
w.measure='normal'; w.absTol=1e-3; w.relTol=1e-3;
q = cubSobol_gCLASS(w);

f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3); exactsol = 1;


% >> q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3); exactsol = 1;
% >> check = abs(exactsol-q) < max(1e-3,1e-3*abs(exactsol))
% check = 1
% 
% 
% Example 3: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [-1,2)^2:
% 
w.f = @(x) exp(-x(:,1).^2-x(:,2).^2); w.domain = [-ones(1,2);2*ones(1,2)];
w.measure='uniform'; w.absTol=1e-3; w.relTol=1e-2;
q = cubSobol_gCLASS(w);
exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;

%
% Example 4: 
% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.
% 
w.f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); w.domain = [-inf(1,1);inf(1,1)];
w.measure='normal'; w.absTol=1e-4; w.relTol=1e-2;
q = cubSobol_gCLASS(w);
price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);

%
% Example 5:
% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the interval
% [0,1)^5 with pure absolute error 1e-5.
% 
w.f = @(x) 8*prod(x,2); w.domain = [zeros(1,5);ones(1,5)];
w.measure='uniform'; w.absTol=1e-5; w.relTol=0;
q = cubSobol_gCLASS(w);
exactsol = 1/4;

%
% Example 6:
% Estimate the integral with integrand f(x) = x1^2+x2^2 over the disk with
% center (0,0) and radius 1 with pure absolute error 1e-4, where x is a vector x = [x1 x2].
% 
w.f = @(x) x(:,1).^2+x(:,2).^2; w.domain = [0,0,1];
w.measure='ball-from-normal'; w.absTol=1e-4; w.relTol=0;
q = cubSobol_gCLASS(w);
exactsol = pi/2;

% Example 7:
% Estimate the integral with integrand f(x) = 10*x1-5*x2^2+x3^3 in the interval [0,2)^3 
% with pure absolute error 1e-6 using two control variates h1(x) = x1 and h2(x) = x2^2.
% 
w.f = @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
w.measure='uniform'; w.absTol=1e-6; w.relTol=0;
w.cv = [8,32/3]; w.domain= [zeros(1,3);2*ones(1,3)];
q = cubSobol_gCLASS(w);
% >> check = abs(exactsol-q) < 1e-6
% check = 1