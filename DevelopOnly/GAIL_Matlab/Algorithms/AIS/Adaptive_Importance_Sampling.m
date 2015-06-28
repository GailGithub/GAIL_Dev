%% Adaptive Importance Sampling (AIS)
%
%
%% Introduction
% Consider the equation below:
%
% $$ \mu = \displaystyle \int_{{R}^d}f(x)dx $$
%
% Where f(x) is a function. If we divide and multiply f(x) by $$ \rho $$(x), we got:
%
% $$ \mu = \displaystyle \int_{{R}^d}\frac{f(x)}{\rho(x)}\rho(x)dx = \int_{{R}^d}g(x)\rho(x)dx $$
%
% Where $$ \rho $$(x) is a probability density function and    
%
% $$ g(x)= \frac{f(x)}{\rho(x)}, x \in R^d $$. 
%
% Importance sampling is the best choice of a probability density function.
%
% The result of the integral 
%
% $$ \int_{{R}^d}g(x)\rho(x)dx $$
%
% Is approximately the mean of the sum of  $$ g(X_i) $$ over the number of samples.
%
% Where X is a random variable defined on $$ R^d $$ with probability
% density function $$ \rho $$(x).
% 
%% Example
%
% Consider that $f(x) = max(x,0)$ and  $$ \rho $$(x) is a normal distribution density function
% given by:
%
% $$ \rho(x)=\frac{1}{\sqrt{2\pi}}exp(\frac{-x^2}{2}) $$
%
% Applying these functions in the integral and simplifying, we obtain: 
%
% $$ integrand = x*exp(\frac{-x^2}{2}) $$
%
% To this integrand we applied a variable transformation $x=t+b$ where 't' is a range of random numbers and 
% 'b' is a determined array (inittialy with three different numbers). From this transformation we obtained:
%
% $$ integrand = (t+b)*exp(-tb-\frac{b^2}{2}) $$
%
% From this integrand we computed the mean (value of the integral)
% using the following programs.


%% meanMC_CLT_AIS
% Our first attempt to use Adaptive Importance Sampling was with the
% program meanMC_CLT_AIS.

type meanMC_CLT_AIS

%% meanMC_CLT_AIS_interp
%
% Next, we developed a new version of the previous program applying a parabolic interpolation to optimize the choice of 'b'
%
type meanMC_CLT_AIS_interp

%% meanMC_AIS_g
%
% Now we are working in the meamMC_g to obtain a guaranteed answer.
%
type meanMC_AIS_g_example

%% Results comparison

disp('For the given function f(x)= max(x,0):')
disp('meanMC_CLT_AIS shows the following results:')
[tmu,out_param]=meanMC_CLT_AIS
disp('meanMC_CLT_AIS_interp shows the following results:')
[tmu,out_param]=meanMC_CLT_AIS_interp
disp('meanMC_AIS_g shows the following results:')
[tmu,out_param]=meanMC_AIS_g_example





