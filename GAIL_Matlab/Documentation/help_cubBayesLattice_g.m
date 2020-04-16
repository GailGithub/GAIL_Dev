%% cubBayesLattice_g
% Bayesian cubature method to estimate the integral
% of a random variable
%% Syntax
% [OBJ,Q] = *cubBayesLattice_g*('f',f,'dim',dim,'absTol',absTol,'relTol',relTol,
% 'order',order,'ptransform',ptransform,'arbMean',arbMean)
%
% [Q,OutP] = *compInteg*(OBJ)
%
%% Description
%
% [OBJ,Q] = *cubBayesLattice_g*('f',f,'dim',d,'absTol',absTol,'relTol',relTol,
% 'order',order,'ptransform',ptransform,'arbMean',arbMean)
% initializes the object with the given parameters and also returns an
% estimate of integral Q.
%
% [Q,OutP] = *compInteg*(OBJ) estimates the integral of f over hyperbox
%   $[0,1]^{d}$ using rank-1 Lattice sampling to within a specified generalized
%   error tolerance, tolfun = max(abstol, reltol*| I |), i.e., | I - Q | <= tolfun
%   with confidence of at least 99%, where I is the true integral value,
%   Q is the estimated integral value, abstol is the absolute error tolerance,
%   and reltol is the relative error tolerance. Usually the reltol determines
%   the accuracy of the estimation, however, if | I | is rather small,
%   then abstol determines the accuracy of the estimation.
%   OutP is the structure holding additional output params, more details provided
%   below. Input f is a function handle that accepts an n x d matrix input,
%   where d is the dimension of the hyperbox, and n is the number of points
%   being evaluated simultaneously.
%
% It is recommended to use COMPINTEG for estimating the integral repeatedly
% after the object initialization.
%
% OutP is the structure holding additional output params, more details
% provided below. Input f is a function handle that accepts an n x d
% matrix input, where d is the dimension of the hyperbox, and n is the
% number of points being evaluated simultaneously.
%
% *Input Arguments*
%
% * f --- the integrand.
%
% * d --- number of dimensions of the integrand.
%
% *Optional Input Arguments*
%
% * absTol --- absolute error tolerance | I - Q | <= absTol. Default is 0.01
%
% * relTol --- relative error tolerance | I - Q | <= I*relTol. Default is 0
%
% * order --- order of the Bernoulli polynomial of the kernel r=1,2.
%             Default is 2
%
% * ptransform --- periodization variable transform to use: 'Baker', 'C0',
%  'C1', 'C1sin', or 'C2sin'. Default is 'C1sin'
%
% * arbMean --- If false, the algorithm assumes the integrand was sampled
%                 from a Gaussian process of zero mean. Default is 'true'
%
% * alpha --- confidence level for a credible interval of Q. Default is 0.01
%
% * mmin --- min number of samples to start with: 2^mmin. Default is 10
%
% * mmax --- max number of samples allowed: 2^mmax. Default is 22
%
% *Output Arguments*
%
% * OutP.n --- number of samples used to compute the integral of f.
%
% * OutP.time --- time to compute the integral in seconds.
%
% <html>
% <ul type="square">
%  <li>OutP.exitFlag --- indicates the exit condition of the
%  algorithm:</li>
%   <ul type="circle">
%                <li>1 - integral computed within the error tolerance and
%                      without exceeding max sample limit 2^mmax </li>
%                <li>2 - used max number of samples and yet not met the
%                      error tolerance</li>
%   </ul>
%  <li>OutP.ErrBd  --- estimated integral error | I - Q |</li>
%  </ul>
% </html>
%
%
%%  Guarantee
%
% This algorithm attempts to calculate the integral of function f over the
% hyperbox $[0,1]^d$ to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level, e.g., 99% when alpha=0.5%. If the
% algorithm terminates without showing any warning messages and provides
% an answer Q, then the following inequality would be satisfied:
%
% Pr(| Q - I | <= tolfun) = 99%
%
% Please refer to our paper for detailed arguments and proofs.
%
%% Examples
%

%%
% *Example 1: Quadratic*
%
% Estimate the integral with integrand $f(x) = x^2$ over the interval $[0,1]$
% with default parameters: order=2, ptransform=C1sin, abstol=0.01, relTol=0

warning('off','GAIL:cubBayesLattice_g:fdnotgiven')
[obj,muhat] = cubBayesLattice_g;
exactInteg = 1.0/3;
warning('on','GAIL:cubBayesLattice_g:fdnotgiven')
check = double(abs(exactInteg-muhat) < 0.01)

%%
% *Example 2: ExpCos*
%
% Estimate the integral with integrand
% $f({x}) = \exp\left(\sum_{i=1}^2cos(2\pi x_i)\right)$ over the
% interval $[0,1]^2$ with parameters: order=2, C1sin variable transform, abstol=0.001,
% relTol=0.01

fun = @(x) exp(sum(cos(2*pi*x), 2));
dim=2; absTol=1e-3; relTol=1e-2;
exactInteg = besseli(0,1)^dim;
inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
inputArgs = [inputArgs {'f',fun, 'dim',dim, 'absTol',absTol,'oneTheta',false}];
obj=cubBayesLattice_g(inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)

%%
% *Example 3: Keister function*
%
% Keister function:
% Estimate the integral with keister function as integrand over the
% interval $[0,1]^2$ with parameters: order=2, C1 variable transform,
% abstol=0.001, relTol=0.01

dim=3; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'f',fKeister,'dim',dim,'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'order',2, 'ptransform','C1','arbMean',true}];
obj=cubBayesLattice_g(inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)

%%
% *Example 4*
%
% Multivariate normal probability:
% Estimate the multivariate normal probability between the hyper interval
% $\left(\begin{array}{c} -6\\ -2\\ -2\end{array}\right) $ and
% $\left(\begin{array}{c} 5\\ 2\\ 1\end{array}\right)$ in $\bf{R}^3$
% having zero mean and covariance
% $\left(\begin{array}{ccc} 4& 1& 1\\ 0& 1& 0.5\\ 0& 0& 0.25 \end{array}\right)$ with
% parameters: order=1, C1sin variable transform,
% abstol=0.001, relTol=0.01

dim=2; absTol=1e-3; relTol=1e-2; fName = 'MVN';
C = [4 1 1; 0 1 0.5; 0 0 0.25]; MVNParams.Cov = C'*C; MVNParams.C = C;
MVNParams.a = [-6 -2 -2]; MVNParams.b = [5 2 1]; MVNParams.mu = 0;
MVNParams.CovProp.C = chol(MVNParams.Cov)';
muBest = 0.676337324357787;
integrand =@(t) GenzFunc(t,MVNParams);
inputArgs={'f',integrand,'dim',dim, 'absTol',absTol,'relTol',relTol};
inputArgs=[inputArgs {'order',1,'ptransform','C1sin','arbMean',true}];
inputArgs=[inputArgs {'useGradient',true}];
[obj,muhat]=cubBayesLattice_g(inputArgs{:});
check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))

%%
% *Example 5: Keister function*
%
% Kernel order r chosen automatically
 
dim=2; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'f',fKeister,'dim',dim,'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'order',0, 'ptransform','C1','arbMean',true}];
obj=cubBayesLattice_g(inputArgs{:});
[muhat,outParams] = compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
check = double(outParams.optParams.r > 0)

%%
% *Example 6*
%
% Another example using dimension specific shape parameter

const = [1E-4 1 1E4];
fun = @(x)sum(bsxfun(@times, const, sin(2*pi*x.^2)), 2);
dim=3; absTol=1e-3; relTol=1e-2;
exactInteg = fresnels(2)*sum(const)/2;
inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
inputArgs = [inputArgs {'f',fun, 'dim',dim, 'absTol',absTol,'oneTheta',false,'useGradient',false}];
obj=cubBayesLattice_g(inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)


%% See Also
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
%% References
%
% [1] Jagadeeswaran Rathinavel, Fred J. Hickernell, Fast automatic Bayesian cubature
%   using lattice sampling.  Stat Comput 29, 1215-1229 (2019).
%   https://doi.org/10.1007/s11222-019-09895-9
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
