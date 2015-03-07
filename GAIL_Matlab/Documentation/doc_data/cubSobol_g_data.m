%%  Guarantee
%
% This algorithm computes the integral of real valued functions in
% dimension d with a prescribed generalized error tolerance. The
% Walsh-Fourier coefficients of the integrand are assumed to be absolutely
% convergent. If the algorithm terminates without warning messages, the
% output is given with guarantees under the assumption that the integrand
% lies inside a cone of functions. The guarantee is based on the decay rate
% of the Walsh-Fourier coefficients. For more details on how the cone is
% defined, please refer to the references below.
%
%% Examples
%
%%
% *Example 1*

% Estimate the integral with integrand f(x) = x1.*x2 in the interval
% [0,1)^2:

  f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,1e-1)

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; hyperbox = [zeros(1,3);ones(1,3)];
  q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3)

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [-1,2)^2:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-3,1e-1)

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [zeros(1,1);ones(1,1)];
  q = cubSobol_g(f,hyperbox,'normal',1e-4,1e-1,'fudge',@(m) 2.^-(2*m))

%%
% *Example 5*

% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the
% interval [0,1)^5 with pure absolute error 1e-5.

  f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0)
