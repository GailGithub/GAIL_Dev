%% Adaptive Importance Sampling (AIS)
%
%% Adaptive Importance Sampling
% Consider the following multivariate integral where _*f(x)*_ is a function:
%
% \begin{align*}
% \mu & =\int_{R^d}^{}{f (x)} \, 
% {\rm
% d}x \\
% \end{align*}
%
% If _*x*_ is a random variable and has a probability density function
% $\rho(x)$, we can assume:
%
% \begin{align*}
% \mu = \int_{{R}^d}\frac{f(x)}{\rho(x)}\rho(x)dx = \int_{{R}^d}g(x)\rho(x)dx 
% \end{align*}
%
% Where:    
%
% \begin{align*}
% g(x)= \frac{f(x)}{\rho(x)}, x \in R^d.
% \end{align*}
%
%
% The integral above can be described then as: 
%
% \begin{align*}
% \mu = \int_{{R}^d}g(x)\rho(x)dx  \rightarrow  \mu\approx
% \frac{1}{n}\sum_{i=1}^{n}{g(X_i)}
% \end{align*}
%
% The importance sampling is the best choice of a probability density function 
% resulting in an alternative distribution whose support is concentrated in the
% truncation region, i.e., sampling more where it matters.

%% Variable Transformation
%
% Suppose that it is complicated to generate _*x*_ with probability density
% $\rho(x)$. We can make use of a variable transformation to obtain a new
% random variable easier to generate and related to _*x*_.
%


%% meanMC_CLT_AIS
%
% Our first algorithm uses the concepts of Adaptive Importance Sample and
% Variable Transformation combined to the Central Limit Theorem to evaluate an integral, searching a value inside
% the transformation that minimizes the variance within an interval.
%
% For example, consider the following integral with a Gaussian probability density $\rho(x)$:
% 
% \begin{align*}
% \mu = \int_{{R}^d}cos(\mid\mid x||)exp(-||x||^2)dx
% \end{align*}
% and
%
% \begin{align*}
% \rho(x) = \frac{exp(-\mid\mid x||^2)/(2b^2)}{((\sqrt{2\pi})b)^d}
% \end{align*}
%
% Calculating the importance sampling formula and using a variable
% transformation _*x = bz*_  where *z* is a random variable and *b* is a fixed value higher than zero, we obtain:
% 
% \begin{align*}
% \mu = \int_{{R}^d}(\sqrt{2\pi}b)^d.cos(b\mid\mid z||)exp((1/2-b^2)\mid\mid z||^2)
% \end{align*}
% 
% Using this function as an input, our algorithm will determine the value
% of *b*, within a determined interval, for which the variance is the
% smallest.
%
% *Example 1:*
d = 3; b = [0.5 2.5]; abstol = 0.002; alpha = 0.01; nSig = 1e4; fudge = 1.2;
Y1=@(z,b) ((sqrt(2.*pi).*b).^d).*cos(b.*sqrt(sum(z.*z,2))).*exp((1/2-b.^2).*sum(z.*z,2));
[tmu,out_param]=meanMC_CLT_AIS(Y1,abstol,alpha,nSig,fudge)

%% meanMC_AIS_g
%
% meanMC_g is a GAIL algorithm that uses Monte Carlo method to estimate the
% mean of a random variable.
% The main input of this algorithm is a function handle that accepts a positive integer input _n_ and
% returns an n x 1 vector of IID instances of the random variable Y.
%
% From our previous algorithm we made the following program, meanMC_AIS_g, using Adaptive Importance Sampling and 
% Variable Transformation to minimize the variance. Once the best value for *b* is found it is used as an input for
% meanMC_g, obtaining a GAIL guaranteed answer.
%
%
%%
% *Program inputs:*
if nargin < 7
   fudge = 1.2; %variance inflation factor
   if nargin < 6
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 5
         alpha = 0.01; %uncertainty
         if nargin < 4
            abstol = 0.002; %absolute error tolerance
            if nargin < 3
                d = 1; %number of dimensions
                if nargin < 2
                    b=[-2 2]; % 'b'interval
                    if nargin < 1 
                        Y1=input('Please inform "g(x)"');% 'g(x)'
                    end                  
                end
            end
         end
      end
   end
end
%%
% *Program restrictions:*
% Checking function input

if isa(Y1,'function_handle') == 0 || nargin(Y1) ~= 2
    warning('meanMC_AIS_g:Y1notafunction',...
    ['"Y" must be a function handle with two variables - "x" and "b".'...
        'A default function "Y(x,b) = randn(n+b,1),b" will be used']);
    Y1 =@(x,b)(randn(x+b,1));
end

% Checking 'b' input

if isa(b,'double') == 0 || numel(b) ~= 2 || issorted(b) == 0
    warning('meanMC_AIS_g:invalidInterval',...
    ['"b" must be an array with 2 elements in ascending order.'...
    'A default interval of [-1 1] will be used.']);
    b=[-1 1];
    
end

% Checking dimension input

if d <= 0 || mod(d,1) ~= 0 
    warning('meanMC_AIS_g:invalidDimension',...
    ['The number of dimensions must be an integer positive number'...
    'A default value d=1 will be used.']);
    d=1;
end

% Checking other inputs

if abstol <= 0 || abstol >1
    warning('meanMC_AIS_g:invalidTolerance',...
        ['The absolut tolerance must be between zero and one'...
    'A default value abstol = 0.002 will be used.']);
    abstol = 0.002;
end

if alpha <= 0
    warning('meanMC_AIS_g:invalidUncertainty',...
        ['The uncertainty must be higher than zero'...
    'A default value alpha = 0.01 will be used.']);
    alpha = 0.01;
end

if nSig <= 0
    warning('meanMC_AIS_g:invalidNumberSamples',...
        ['The number of samples must be higher than zero'...
    'A default value nSig = 1e4 will be used.']);
    nSig = 1e4;
end

if fudge <= 0
    warning('meanMC_AIS_g:invalidInflationFactor',...
        ['The inflation factor must be higher than zero'...
    'A default value fudge = 1.2 will be used.']);
    fudge = 0.01;
end

%%
% *Main Structure:*

out_param_AIS.alpha = alpha; % Save the input parameters to a structure.
out_param_AIS.abstol = abstol;
out_param_AIS.fudge = fudge;
out_param_AIS.nSig = nSig;

b_vec=[b(1),((b(1)+b(2))/2),b(2)]; % Generates a vector with 3 values equally spaced
%within the interval defined.

Y = @(n,b)Y1(randn(n,d),b); % Integrand evaluated at the sample points. 

tstart=tic; % Starts the clock.

var_b=b_vec; % Avoids the change of vector size inside the loop.

% Checking the variance for each element in b_vec:
for i=1:numel(b_vec)
    b_value=b_vec(i);    
    var_b(i)=var(Y(nSig,b_value));
end
[S_var,S_pos]=min(var_b); % Computes the position and the value of the smaller
% variance calculated from 'b_vec'

% Parabolic interpolation between b_vec and calculated variance:
A=[b_vec'.^2 b_vec' ones(3,1)];
p=A\var_b';

% Minimum search using the approximated parabola:
[x]=(-p(2)/(2*p(1)));
% Variance calculation using the value estimated as the minimum 'x'
var_bx=var(Y(nSig,x)); 


% Checking the best value for b:
if var_bx < S_var && var_bx > 0
    out_param_AIS.b_value = x;
    out_param_AIS.var = var_bx;
else
    out_param_AIS.b_value = b_vec(S_pos);
    out_param_AIS.var = S_var;
end
out_param_AIS.time=toc(tstart); %elapsed time
%%
% *meanMC_g calculation:*

[tmu, out_param_MCg]=meanMC_g(@(n)Y(n,out_param_AIS.b_value),out_param_AIS.abstol,0,out_param_AIS.alpha,out_param_AIS.fudge,out_param_AIS.nSig);

out_param_AIS.nTotal= 4.*nSig+(out_param_MCg.ntot);%total number of samples used
out_param_AIS.sig0 = sqrt(out_param_AIS.var); %standard deviation

%% 
% *Example 2:* Using the same inputs of the Example 1:

[tmu,out_param_AIS, out_param_MCg]=meanMC_AIS_g(Y1,b,d,abstol,alpha,nSig,fudge)

%% References
%
% 
% * HICKERNELL, F. J. _*Monte Carlo and Quasi-Monte Carlo Methods*_. Illinois
% Institute of Technology May, 2015.
% 
%
%%
%
% *Authors*
%
% BRITO, Rafael de Miranda.
% DE PAULO, Ricardo Freitas.
% SABARENSE, Mariane de Carvalho.
