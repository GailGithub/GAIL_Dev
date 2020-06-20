%% cubBayesLattice_g
% Bayesian cubature method to estimate the integral of a random variable 
% using rank-1 Lattices over a d-dimensional region within a 
% specified generalized error tolerance with guarantees under Bayesian
% assumptions.
%% Syntax
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,'absTol',absTol,'relTol',relTol,
%           'order',order,'ptransform',ptransform,'arbMean',arbMean)
%
% [OBJ] = *cubBayesLattice_g*(f,dim,'absTol',absTol,'relTol',relTol,
%         'order',order,'ptransform',ptransform,'arbMean',arbMean)
%
% [Q,OutP] = *compInteg*(OBJ)
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim)
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,absTol,relTol)
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,inParams)
%
%% Description
% 
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,'absTol',absTol,'relTol',relTol,
%   'order',order,'ptransform',ptransform,'arbMean',arbMean) Initializes
%   the object with the given parameters and also returns an
%   estimate of integral Q.
%
% [Q,OutP] = *compInteg*(OBJ) estimates the integral of f over hyperbox
%   $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling to within a specified generalized
%   error tolerance, tolfun = max(abstol, reltol*| I |), i.e., | I - Q | <=
%   tolfun with confidence of at least 99%, where I is the true integral
%   value, Q is the estimated integral value, abstol is the absolute error
%   tolerance, and reltol is the relative error tolerance. Usually the
%   reltol determines the accuracy of the estimation; however, if | I | is
%   rather small, then abstol determines the accuracy of the estimation.
%   Given the construction of our Lattices, d must be a positive integer
%   with 1 <= dim <= 600. For higher dimensions, it is recommended to use 
%   simpler periodization transformation like 'Baker'.
%
% It is recommended to use *compInteg* for estimating the integral repeatedly
% after the object initialization.
% 
% OutP is the structure holding additional output params, more details
% provided below. Input f is a function handle that accepts an n x d
% matrix input, where d is the dimension of the hyperbox, and n is the
% number of points being evaluated simultaneously.
%
% The following additional input parameter passing styles also supported:
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim) estimates the integral of f over
%   hyperbox $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling. All other input parameters
%   are initialized with default values as given below. Returns the initialized
%   object OBJ and the estimate of integral Q.
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,absTol,relTol); estimates the integral 
%   of f over hyperbox $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling. All parameters
%   should be input in the order specified above. The answer is given within 
%   the generalized error tolerance tolfun. All other input parameters
%    are initialized with default values as given below.  
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,inParms); estimates the integral 
%   of f over hyperbox $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling. 
%   The structure inParams shall hold the optional input parameters.
%
% *Input Arguments*
%
% * f --- the integrand.
%
% * dim --- number of dimensions of the integrand.
%
% *Optional Input Arguments*
%
% * absTol --- absolute error tolerance | I - Q | <= absTol. Default is 0.01
%
% * relTol --- relative error tolerance | I - Q | <= I*relTol. Default is 0
%
% * order --- order of the Bernoulli polynomial of the kernel r=1,2. 
%             If r==0, algorithm automatically chooses the kernel order
%             which can be a non-integer value.
%             Default is 2
%
% * ptransform --- periodization variable transform to use: 'Baker', 'C0',
%                  'C1', 'C1sin', or 'C2sin'. Default is 'C1sin'
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
% <html>
% <ul type="square">
% <li>stopCriterion --- stopping criterion to use. Supports three options: </li>
%   <ul type="circle">
%                <li>1) MLE: Empirical Bayes</li>
%                <li>2) GCV: Generalized Cross Validation</li>
%                <li>3) full: Full Bayes</li>
%   </ul>
%    Default is MLE: Empirical Bayes
%  </ul>
% </html>
%
% * useGradient --- If true uses gradient descent in parameter search.
%                   Default is false
%
% * oneTheta --- If true uses common shape parameter for all dimensions,
%                 else allow shape parameter vary across dimensions.
%                 Default is true
%
% *Output Arguments*
%
% * n --- number of samples used to compute the integral of f.
%
% * time --- time to compute the integral in seconds.
%
% <html>
% <ul type="square">
%  <li>exitFlag --- indicates the exit condition of the
%  algorithm:</li>
%   <ul type="circle">
%                <li>1 - integral computed within the error tolerance and
%                      without exceeding max sample limit 2^mmax </li>
%                <li>2 - used max number of samples and yet not met the
%                      error tolerance</li>
%   </ul>
%  </ul>
% </html>
%
% * ErrBd  --- estimated integral error | I - Q |
%
% * optParams --- optional parameters useful to debug and get better
%                  understanding of the algorithm
%
% * optParams.aMLEAll --- returns the shape parameters computed
%
%
%
%%  Guarantee
%
% This algorithm attempts to calculate the integral of function f over the
% hyperbox $[0,1]^{\mathrm{dim}}$ to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level, e.g., 99% when alpha=0.5%. If the
% algorithm terminates without showing any warning messages and provides
% an answer Q, then the following inequality would be satisfied:
%
% Pr(| Q - I | <= tolfun) = 99%
%
% Please refer to our paper [1] for detailed arguments and proofs.
%
%% Examples
%

%%
% *Example 1: Integrating a simple Quadratic function*
%
% Estimate the integral with integrand $f(x) = x^2$ over the interval $[0,1]$
% with default parameters: order=2, ptransform=C1sin, abstol=0.01, relTol=0

warning('off','GAIL:cubBayesLattice_g:fdnotgiven')
[~,muhat] = cubBayesLattice_g;
exactInteg = 1.0/3;
warning('on','GAIL:cubBayesLattice_g:fdnotgiven')
check = double(abs(exactInteg-muhat) < 0.01)

%%
% *Example 2: ExpCos*
%
% Estimate the integral of Exponential of Cosine function
% $f({x}) = \exp\left(\sum_{i=1}^2cos(2\pi x_i)\right)$ over the
% interval $[0,1]^2$ with parameters: order=2, C1sin variable transform, abstol=0.001,
% relTol=0.01

fun = @(x) exp(sum(cos(2*pi*x), 2));
dim=2; absTol=1e-3; relTol=1e-2;
exactInteg = besseli(0,1)^dim;
inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
inputArgs = [inputArgs {'absTol',absTol,'oneTheta',false}];
obj=cubBayesLattice_g(fun,dim,inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)

%%
% *Example 3: Keister function*
%
% Estimate the integral with keister function as integrand over the
% interval $[0,1]^2$ with parameters: order=2, C1 variable transform,
% abstol=0.001, relTol=0.01

dim=3; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'order',2, 'ptransform','C1','arbMean',true}];
obj=cubBayesLattice_g(fKeister,dim,inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)

%%
% *Example 4:  Multivariate normal probability*
% 
% Estimate the multivariate normal probability for the given hyper interval
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
inputArgs={'absTol',absTol,'relTol',relTol};
inputArgs=[inputArgs {'order',1,'ptransform','C1sin','arbMean',true}];
inputArgs=[inputArgs {'useGradient',true}];
[~,muhat]=cubBayesLattice_g(integrand,dim, inputArgs{:});
check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))

%%
% *Example 5: Keister function*
%
% Estimating the Keister integral with Kernel order r chosen automatically
 
dim=2; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'order',0, 'ptransform','C1','arbMean',true}];
obj=cubBayesLattice_g(fKeister,dim,inputArgs{:});
[muhat,outParams] = compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
check = double(outParams.optParams.r > 0)

%%
% *Example 6*
%
% A simple example which uses dimension specific shape parameter

const = [1E-4 1 1E4];
fun = @(x)sum(bsxfun(@times, const, sin(2*pi*x.^2)), 2);
dim=3; absTol=1e-3; relTol=1e-2;
exactInteg = fresnels(2)*sum(const)/2;
inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
inputArgs = [inputArgs {'absTol',absTol,'oneTheta',false,'useGradient',false}];
obj=cubBayesLattice_g(fun, dim, inputArgs{:});
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
% <html>
% <a href="help_cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Jagadeeswaran Rathinavel, Fred J. Hickernell, Fast automatic Bayesian 
%   cubature using lattice sampling.  Stat Comput 29, 1215-1229 (2019).
%   https://doi.org/10.1007/s11222-019-09895-9
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using
%   matching kernels and designs," PhD thesis, Illinois Institute of Technology, 2019.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
