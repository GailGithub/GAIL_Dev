%%  Guarantee
% This algorithm attempts to calculate the integral of function f over a
% hyperbox to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level 1-alpha. If the algorithm terminated
% without showing any warning messages and provide an answer Q, then the
% follow inequality would be satisfied:
% 
% Pr(| Q - I | <= tolfun) >= 1-alpha
%
% The cost of the algorithm, N_tot, is also bounded above by N_up, which is
% a function in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta.
% And the following inequality holds:
% 
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
%
%%  Examples
% *Example 1*

% Estimate the integral with integrand f(x) = sin(x) over the interval
% [1;2]
% 

 f=@(x) sin(x);interval = [1;2];
 Q = cubMC_g(f,interval,'uniform',1e-3,1e-2)
 
%% 
% *Example 2*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) over the
% hyperbox [0 0;1 1], where x is a vector x = [x1 x2].
% 

 f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [0 0;1 1];
 Q = cubMC_g(f,hyperbox,'measure','uniform','abstol',1e-3,...
     'reltol',1e-13)

%%
% *Example 3*

% Estimate the integral with integrand f(x) = 2^d*prod(x1*x2*...*xd) +
% 0.555 over the hyperbox [zeros(1,d);ones(1,d)], where x is a vector x =
% [x1 x2... xd].


  d=3;f=@(x) 2^d*prod(x,2)+0.555; hyperbox =[zeros(1,d);ones(1,d)];
  in_param.abstol = 1e-3;in_param.reltol=1e-3;
  Q = cubMC_g(f,hyperbox,in_param)

%%
% *Example 4* 

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [-inf -inf;inf inf], where x is a vector x = [x1 x2].


 f=@(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-inf -inf;inf inf];
 Q = cubMC_g(f,hyperbox,'normal',0,1e-2)
