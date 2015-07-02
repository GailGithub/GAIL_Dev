%% Adaptive Importance Sampling (AIS)
%
%% Adaptive Importance Sampling
% Consider the following multivariate integral,
%
% \begin{equation*}
% \mu =\int_{\mathbb{R}^d}^{}{f (\boldsymbol{x})} \, 
% {\rm
% d}\boldsymbol{x} \\
% \end{equation*}
%
% where \(f\) is some function defined on the \(d\)-dimensional real
% vectors.  If \(\boldsymbol{X}\) is a random variable with probability density
% function $\rho$, then $\mu$ can be represented as the mean of a function
% of \(\boldsymbol{X}\).
%
% \begin{equation*} \mu =
% \int_{\mathbb{R}^d}\frac{f(\boldsymbol{x})}{\rho(\boldsymbol{x})}
% \rho(\boldsymbol{x}) \, {\rm d} \boldsymbol{x} =
% \int_{\mathbb{R}^d}g(\boldsymbol{x}) \rho(\boldsymbol{x}) \, {\rm d}
% \boldsymbol{x} =: \mathbb{E}[g(\boldsymbol{X})], \end{equation*}
%
% where   
%
% \begin{equation*} g(\boldsymbol{x}) := \frac{f(\boldsymbol{x})}
% {\rho(\boldsymbol{x})}, \qquad \boldsymbol{x} \in R^d. \end{equation*}
%
% For independent and identically distributed (IID) \(\boldsymbol{X}_i\)
% the integral may be approximated by a sum: 
%
% \begin{equation*} \lim_{n \to \infty} \mathbb{E}[\lvert \mu - \hat{\mu}_n
% \rvert^2] = 0, \qquad \mbox{where } \hat{\mu}_n :=
% \frac{1}{n}\sum_{i=1}^{n} g(X_i)  
% \end{equation*}
%
% Note that this is true for all choices of \(\rho\).  Importance sampling
% chooses \(\rho\) to make \(\mathbb{E}[\lvert \mu - \hat{\mu}_n
% \rvert^2]\) tend to zero faster as \(n \to \infty\), i.e., sampling more
% where it matters.

%% Variable Transformation
%
% Suppose that it is complicated to generate \(\boldsymbol{X}\) with
% probability density \(\rho\). We can make use of a variable
% transformation to obtain a new random variable easier to generate and
% related to \(\boldsymbol{X}\).
%


%% meanMC_CLT_AIS
%
% Our first algorithm uses the concepts of Adaptive Importance Sample and
% Variable Transformation combined with the Central Limit Theorem to
% evaluate an integral, searching a value inside the transformation that
% minimizes the variance within an interval.
%
% For example, consider the following integral and Gaussian probability
% density $\rho$:
% 
% \begin{equation*} 
% \mu = \int_{\mathbb{R}^d} \cos(\lVert \boldsymbol{x} \rVert)
% \exp(-\lVert \boldsymbol{x} \rVert^2) \, {\rm d} \boldsymbol{x} \qquad \mbox{and} \qquad
% \rho(\boldsymbol{x}) = \frac{\exp(-\lVert \boldsymbol{x}
% \rVert^2)/(2b^2)}{(\sqrt{2\pi}b)^d} 
% \end{equation*}
%
% Calculating the importance sampling formula and using a variable
% transformation \(\boldsymbol{x} = b\boldsymbol{z}\)  where \(b>0\), we
% obtain:
% 
% \begin{equation*} \mu = \int_{\mathbb{R}^d}(\sqrt{2\pi}b)^d \cos(b\lVert
% \boldsymbol{z} \rVert) \exp((1/2-b^2)\lVert \boldsymbol{z} \rVert^2) \,
% \frac{\exp(-\lVert \boldsymbol{z} \rVert^2)/2}{(\sqrt{2\pi})^d} \,
% {\rm d} \boldsymbol{z} \end{equation*}
% 
% Using this function as an input, our algorithm will determine the value
% of *b*, within a determined interval, for which the variance is the
% smallest.

%% Example 1:
b = [0.5 2.5]; abstol = 0.002; alpha = 0.01; nSig = 1e4; fudge = 1.2; d=1;
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

type meanMC_AIS_g

%% Example 2: 
% Using the same function of the Example 1:

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
