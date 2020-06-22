function [q,out_param,y,kappanumap] = cubLattice_g(varargin)
%CUBLATTICE_G Quasi-Monte Carlo method using rank-1 Lattices cubature over
%a d-dimensional region to integrate within a specified generalized error
%tolerance with guarantees under Fourier coefficients cone decay
%assumptions.
%
%   [q,out_param] = CUBLATTICE_G(f,hyperbox) estimates the integral of f
%   over the d-dimensional region described by hyperbox, and with an error
%   guaranteed not to be greater than a specific generalized error
%   tolerance, tolfun:=max(abstol,reltol*| integral(f) |). Input f is a
%   function handle. f should accept an n x d matrix input, where d is the
%   dimension and n is the number of points being evaluated simultaneously.
%
%   When measure is 'uniform', the input hyperbox is a 2 x d matrix, where
%   the first row corresponds to the lower limits and the second row
%   corresponds to the upper limits of the integral. When measure is
%   'uniform ball' or 'uniform sphere', the input hyperbox is a vector with
%   d+1 elements, where the first d values correspond to the center of the
%   ball and the last value corresponds to the radius of the ball. For
%   these last two measures, a user can optionally specify what
%   transformation should be used in order to get a uniform distribution on
%   a ball. When measure is 'uniform ball_box', the box-to-ball
%   transformation, which gets a set of points uniformly distributed on a
%   ball from a set of points uniformly distributed on a box, will be used.
%   When measure is 'uniform ball_normal', the normal-to-ball
%   transformation, which gets a set of points uniformly distributed on a
%   ball from a set of points normally distributed on the space, will be
%   used. Similarly, the measures 'uniform sphere_box' and 'uniform
%   sphere_normal' can be used to specify the desired transformations. The
%   default transformations are the box-to-ball and the box-to-sphere
%   transformations, depending on the region of integration. Given the
%   construction of our Lattices, d must be a positive integer with 1 <= d
%   <= 600.
%
%   q = CUBLATTICE_G(f,hyperbox,measure,abstol,reltol)
%   estimates the integral of f over the hyperbox. The answer
%   is given within the generalized error tolerance tolfun. All parameters
%   should be input in the order specified above. If an input is not specified,
%   the default value is used. Note that if an input is not specified,
%   the remaining tail cannot be specified either. Inputs f and hyperbox
%   are required. The other optional inputs are in the correct order:
%   measure,abstol,reltol,shift,mmin,mmax,fudge, and transform.
%
%   q = CUBLATTICE_G(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%   estimates the integral of f over the hyperbox. The answer
%   is given within the generalized error tolerance tolfun. All the field-value
%   pairs are optional and can be supplied in any order. If an input is not
%   specified, the default value is used.
%
%   q = CUBLATTICE_G(f,hyperbox,in_param) estimates the integral of f over the
%   hyperbox. The answer is given within the generalized error tolerance tolfun.
%
%   Input Arguments
%
%     f --- the integrand whose input should be a matrix n x d where n is
%     the number of data points and d the dimension, which cannot be
%     greater than 600. By default f is f=@ x.^2.
%
%     hyperbox --- the integration region defined by its bounds. When
%     measure is 'uniform' or 'normal', hyperbox must be a 2 x d matrix,
%     where the first row corresponds to the lower limits and the second
%     row corresponds to the upper limits of the integral. When measure is
%     'uniform ball' or 'uniform sphere', the input hyperbox is a vector
%     with d+1 elements, where the first d values correspond to the center
%     of the ball and the last value corresponds to the radius of the ball.
%     The default value is [0;1].
%
%     in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%     measure of a uniformly distributed random variable in the hyperbox or
%     normally distributed with covariance matrix I_d. The possible values
%     are 'uniform', 'normal', 'uniform ball', 'uniform ball_box', 'uniform
%     ball_normal', 'uniform sphere', 'uniform sphere_box' and 'uniform
%     sphere_normal'. For 'uniform', the hyperbox must be a finite volume,
%     for 'normal', the hyperbox can only be defined as (-Inf,Inf)^d and,
%     for 'uniform ball' or 'uniform sphere', hyperbox must have finite
%     values for the coordinates of the center and a finite positive value
%     for the radius. By default it is 'uniform'.
%
%     in_param.abstol --- the absolute error tolerance, abstol>=0. By
%     default it is 1e-4. For pure absolute tolerance, set in_param.reltol
%     = 0.
%
%     in_param.reltol --- the relative error tolerance, which should be in
%     [0,1]. Default value is 1e-2. For pure absolute tolerance, set
%     in_param.abstol = 0.
%
%   Optional Input Arguments
%
%     in_param.shift --- the Rank-1 lattices can be shifted to avoid the
%     origin or other particular points. The shift is a vector in [0,1]^d.
%     By default we consider a shift uniformly sampled from [0,1]^d.
%
%     in_param.mmin --- the minimum number of points to start is 2^mmin.
%     The cone condition on the Fourier coefficients decay requires a
%     minimum number of points to start. The advice is to consider at least
%     mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%     default it is 10.
%
%     in_param.mmax --- the maximum budget is 2^mmax. By construction of
%     our Lattices generator, mmax is a positive integer such that
%     mmin<=mmax. mmax should not be bigger than the gail.lattice_gen
%     allows. The default value is 20.
%
%     in_param.fudge --- the positive function multiplying the finite
%     sum of Fast Fourier coefficients specified in the cone of functions.
%     This input is a function handle. The fudge should accept an array of
%     nonnegative integers being evaluated simultaneously. For more
%     technical information about this parameter, refer to the references.
%     By default it is @(m) 5*2.^-m.
%
%     in_param.transform --- the algorithm is defined for continuous
%     periodic functions. If the input function f is not, there are 5
%     types of transform to periodize it without modifying the result.
%     By default it is the Baker's transform. The options are:
%       id : no transformation.
%       Baker : Baker's transform or tent map in each coordinate. Preserving
%                 only continuity but simple to compute. Chosen by default.
%       C0 : polynomial transformation only preserving continuity.
%       C1 : polynomial transformation preserving the first derivative.
%       C1sin : Sidi's transform with sine, preserving the first derivative.
%                 This is in general a better option than 'C1'.
%
%   Output Arguments
%
%     q --- the estimated value of the integral.
%
%     out_param.d --- dimension over which the algorithm integrated.
%
%     out_param.n --- number of Rank-1 lattice points used for computing
%     the integral of f.
%
%     out_param.bound_err --- predicted bound on the error based on the
%     cone condition. If the function lies in the cone, the real error will
%     be smaller than generalized tolerance.
%
%     out_param.time --- time elapsed in seconds when calling cubLattice_g.
%
%     out_param.exitflag --- this is a binary vector stating whether
%     warning flags arise. These flags tell about which conditions make the
%     final result certainly not guaranteed. One flag is considered arisen
%     when its value is 1. The following list explains the flags in the
%     respective vector order:
%
%                       1 : If reached overbudget, meaning the max budget
%                       is attained without reaching the guaranteed error
%                       tolerance.
%
%                       2 : If the function lies outside the cone, results
%                       are not guaranteed to be accurate. Note that this
%                       parameter is computed on the transformed function,
%                       not the input function. For more information on the
%                       transforms, check the input parameter
%                       in_param.transform; for information about the cone
%                       definition, check the article mentioned below.
%
%  Guarantee
% This algorithm computes the integral of real valued functions in [0,1]^d
% with a prescribed generalized error tolerance. The Fourier coefficients
% of the integrand are assumed to be absolutely convergent. If the
% algorithm terminates without warning messages, the output is given with
% guarantees under the assumption that the integrand lies inside a cone of
% functions. The guarantee is based on the decay rate of the Fourier
% coefficients. For integration over domains other than [0,1]^d, this cone
% condition applies to f \circ \psi (the composition of the
% functions) where \psi is the transformation function for [0,1]^d to
% the desired region. For more details on how the cone is defined, please
% refer to the references below.
%
%  Examples
%
% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1]^2:
%
% >> f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)];
% >> q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','C1sin'); exactsol = 1/4;
% >> check = double(abs(exactsol-q) < 1e-5)
% check = 1
%
%
% Example 2:
% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2 in the
% interval R^3 where x1, x2 and x3 are normally distributed:
%
% >> f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
% >> q = cubLattice_g(f,hyperbox,'normal',1e-3,1e-3,'transform','C1sin','shift',2^(-25)*ones(1,3)); exactsol = 1;
% >> check = double(abs(exactsol-q) < max(1e-3,1e-3*abs(exactsol)))
% check = 1
%
%
% Example 3:
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [-1,2]^2:
%
% >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
% >> q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1'); exactsol = 1/9*(sqrt(pi)/2*(erf(2)+erf(1)))^2;
% >> check = double(abs(exactsol-q) < max(1e-3,1e-2*abs(exactsol)))
% check = 1
%
%
% Example 4:
% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.
%
% >> f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [-inf(1,1);inf(1,1)];
% >> q = cubLattice_g(f,hyperbox,'normal',1e-4,1e-2,'transform','C1sin'); price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
% >> check = double(abs(price-q) < max(1e-4,1e-2*abs(price)))
% check = 1
%
%
% Example 5:
% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the
% interval [0,1]^5 with pure absolute error 1e-5.
%
% >> f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
% >> q = cubLattice_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
% >> check = double(abs(exactsol-q) < 1e-5)
% check = 1
%
%
% Example 6:
% Estimate the integral with integrand f(x) = 3./(5-4*(cos(2*pi*x))) in the
% interval [0,1] with pure absolute error 1e-5.
%
% >> f = @(x) 3./(5-4*(cos(2*pi*x))); hyperbox = [0;1];
% >> q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','id'); exactsol = 1;
% >> check = double(abs(exactsol-q) < 1e-5)
% check = 1
%
%
% Example 7:
% Estimate the integral with integrand f(x) = x1^2+x2^2 over the disk with
% center (0,0) and radius 1 with pure absolute error 1e-4, where x is a
% vector x = [x1 x2].
%
% >> f = @(x) x(:,1).^2+x(:,2).^2; hyperbox = [0,0,1];
% >> q = cubLattice_g(f,hyperbox,'uniform ball','abstol',1e-4,'reltol',0); exactsol = pi/2;
% >> check = double(abs(exactsol-q) < 1e-4)
% check = 1
%
%
%   See also CUBSOBOL_G, CUBMC_G, MEANMC_G, INTEGRAL_G, CUBBAYESLATTICE_G
%
%  References
%
%   [1] Lluis Antoni Jimenez Rugama and Fred J. Hickernell, "Adaptive
%   multidimensional integration based on rank-1 lattices," Monte Carlo and
%   Quasi-Monte Carlo  Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools
%   and D. Nuyens, eds.), Springer Proceedings in Mathematics and Statistics,
%   vol. 163, Springer-Verlag, Berlin, 2016, arXiv:1411.1966, pp. 407-422.
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice and Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   [6] Kai-Tai Fang and Yuan Wang, Number-theoretic Methods in 
% Statistics, Chapman & Hall, London, 1994.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

t_start = tic;
%% Initial important cone factors and Check-initialize parameters
r_lag = 4; %distance between coefficients summed and those computed
[f,hyperbox,out_param] = cubLattice_g_param(r_lag,varargin{:});

%------------------------------------------------------------------------------
% TRANSFORMATION

%changing the integrand and the hyperbox when measure is uniform ball or
%sphere by applying the appropriate transformation
if strcmpi(out_param.measure,'uniform ball') || strcmpi(out_param.measure,'uniform sphere')% using uniformly distributed samples on a ball or sphere
    if strcmp(out_param.measure,'uniform sphere') && out_param.transf == 1 %box-to-sphere transformation
        out_param.d = out_param.d + 1; % changing out_param.d to the dimension of the sphere
        out_param.shift = [out_param.shift rand];
    end

    if strcmpi(out_param.measure,'uniform ball')% using the formula of the volume of a ball
        volume = ((2.0*pi^(out_param.d/2.0))/(out_param.d*gamma(out_param.d/2.0)))*out_param.radius^out_param.d; %volume of a d-dimentional ball
    else % using the formula of the volume of a sphere
        volume = ((2.0*pi^(out_param.d/2.0))/(gamma(out_param.d/2.0)))*out_param.radius^(out_param.d - 1); %volume of a d-dimentional sphere
    end

    if out_param.transf == 1 % box-to-ball or box-to-sphere transformation should be used
        if out_param.d == 1 % It is not necessary to multiply the function f by the volume, since no transformation is being made
            hyperbox = [hyperbox - out_param.radius; hyperbox + out_param.radius];% for one dimension, the ball is actually an interval
            out_param.measure = 'uniform';% then a uniform distribution on a box can be used
        else
            if strcmpi(out_param.measure,'uniform ball') % box-to-ball transformation
                f = @(t) f(gail.domain_balls_spheres.ball_psi_1(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
            else %  % box-to-sphere transformation
                f = @(t) f(gail.domain_balls_spheres.sphere_psi_1(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
                out_param.d = out_param.d - 1;% the box-to-sphere transformation takes points from a (d-1)-dimensional box to a d-dimensional sphere
                out_param.shift = out_param.shift(1:end-1);
            end
            hyperbox = [zeros(1, out_param.d); ones(1, out_param.d)];% the hyperbox must be the domain of the transformation, which is a unit box
            out_param.measure = 'uniform';% then a uniform distribution on a box can be used
        end
    else % normal-to-ball or normal-to-sphere transformation should be used
        if strcmpi(out_param.measure,'uniform ball') % normal-to-ball transformation
            f = @(t) f(gail.domain_balls_spheres.ball_psi_2(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
        else % normal-to-sphere transformation
            f = @(t) f(gail.domain_balls_spheres.sphere_psi_2(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
        end
        hyperbox = bsxfun(@plus, zeros(2, out_param.d), [-inf; inf]);% the hyperbox must be the domain of the transformation, which is a this unit box
        out_param.measure = 'normal';% then a normal distribution can be used
    end
end

%------------------------------------------------------------------------------

l_star = out_param.mmin - r_lag; % Minimum gathering of points for the sums of DFT
omg_circ = @(m) 2.^(-m);
omg_hat = @(m) out_param.fudge(m)/((1+out_param.fudge(r_lag))*omg_circ(r_lag));

if strcmp(out_param.measure,'normal')
   f=@(x) f(gail.stdnorminv(x));
elseif strcmp(out_param.measure,'uniform')
   f=@(x) f(bsxfun(@plus,hyperbox(1,:),bsxfun(@times,(hyperbox(2,:)-hyperbox(1,:)),x))); % a + (b-a)x = u
end

if strcmp(out_param.transform,'Baker')
    f=@(x) f(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(out_param.transform,'C0')
    f=@(x) f(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(out_param.transform,'C1')
    f=@(x) f(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(out_param.transform,'C1sin')
    f=@(x) f(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
end

%% Main algorithm - Preallocation
Stilde=zeros(out_param.mmax-out_param.mmin+1,1); %initialize sum of DFT terms
CStilde_low = -inf(1,out_param.mmax-l_star+1); %initialize various sums of DFT terms for necessary conditions
CStilde_up = inf(1,out_param.mmax-l_star+1); %initialize various sums of DFT terms for necessary conditions
errest=zeros(out_param.mmax-out_param.mmin+1,1); %initialize error estimates
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1); %initialize approximations to integral
exit_len = 2;
out_param.exit=false(1,exit_len); %we start the algorithm with all warning flags down

%% Initial points and FFT
out_param.n=2^out_param.mmin; %total number of points to start with
n0=out_param.n; %initial number of points
xpts=mod(bsxfun(@plus, gail.lattice_gen(1,n0,out_param.d), out_param.shift),1); %grab Lattice points
y=f(xpts); %evaluate integrand
yval=y;

%% Compute initial FFT
for l=0:out_param.mmin-1
   nl=2^l;
   nmminlm1=2^(out_param.mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;
end
% y now contains the FFT coefficients

%% Approximate integral
q=mean(yval);
appxinteg(1)=q;

%% Create kappanumap implicitly from the data
kappanumap=(1:out_param.n)'; %initialize map
for l=out_param.mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
   flip=find(newone>oldone); %which in the pair are the larger ones
   if ~isempty(flip)
       flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
       flipall=flipall(:);
       temp=kappanumap(nl+1+flipall); %then flip
       kappanumap(nl+1+flipall)=kappanumap(1+flipall); %them
       kappanumap(1+flipall)=temp; %around
   end
end

%% Compute Stilde
nllstart=int64(2^(out_param.mmin-r_lag-1));
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
out_param.bound_err=out_param.fudge(out_param.mmin)*Stilde(1);
errest(1)=out_param.bound_err;

% Necessary conditions
for l = l_star:out_param.mmin % Storing the information for the necessary conditions
    C_low = 1/(1+omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    C_up = 1/(1-omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    CStilde_low(l-l_star+1) = max(CStilde_low(l-l_star+1),C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
    if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
        CStilde_up(l-l_star+1) = min(CStilde_up(l-l_star+1),C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
    end
end
if any(CStilde_low(:) > CStilde_up(:))
   out_param.exit(2) = true;
end

% Check the end of the algorithm
q = q - errest(1)*(max(out_param.abstol, out_param.reltol*abs(q + errest(1)))...
    - max(out_param.abstol, out_param.reltol*abs(q - errest(1))))/...
    (max(out_param.abstol, out_param.reltol*abs(q + errest(1)))...
    + max(out_param.abstol, out_param.reltol*abs(q - errest(1)))); % Optimal estimator
appxinteg(1)=q;

is_done = false;
if 4*errest(1)^2/(max(out_param.abstol, out_param.reltol*abs(q + errest(1)))...
    + max(out_param.abstol, out_param.reltol*abs(q - errest(1))))^2 <= 1
   out_param.time=toc(t_start);
   is_done = true;
elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
   out_param.exit(1) = true;
   is_done = true;
end

%% Loop over m
for m=out_param.mmin+1:out_param.mmax
   if is_done,
       break;
   end
   out_param.n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=mod(bsxfun(@plus, gail.lattice_gen(nnext+1,2*nnext,out_param.d), out_param.shift),1);
   n0=n0+nnext;
   ynext=f(xnext);
   yval=[yval; ynext]; %#ok<*AGROW>

   %% Compute initial FFT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
      coefv=repmat(coef,nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+coefv.*oddval)/2;
      ynext(~ptind)=(evenval-coefv.*oddval)/2;
   end

   %% Compute FFT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;

   %% Update kappanumap
   kappanumap=[kappanumap; 2^(m-1)+kappanumap]; %initialize map
   for l=m-1:-1:m-r_lag
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
      flip=find(newone>oldone);
      if ~isempty(flip)
          flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
          flipall=flipall(:);
          temp=kappanumap(nl+1+flipall);
          kappanumap(nl+1+flipall)=kappanumap(1+flipall);
          kappanumap(1+flipall)=temp;
      end
   end

   %% Compute Stilde
   nllstart=int64(2^(m-r_lag-1));
   meff=m-out_param.mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   out_param.bound_err=out_param.fudge(m)*Stilde(meff);
   errest(meff)=out_param.bound_err;

   % Necessary conditions
   for l = l_star:m % Storing the information for the necessary conditions
        C_low = 1/(1+omg_hat(m-l)*omg_circ(m-l));
        C_up = 1/(1-omg_hat(m-l)*omg_circ(m-l));
        CStilde_low(l-l_star+1) = max(CStilde_low(l-l_star+1),C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
        if (omg_hat(m-l)*omg_circ(m-l) < 1)
            CStilde_up(l-l_star+1) = min(CStilde_up(l-l_star+1),C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
        end
   end

   if any(CStilde_low(:) > CStilde_up(:))
       out_param.exit(2) = true;
   end

   %% Approximate integral
   q=mean(yval);

   % Check the end of the algorithm
   q = q - errest(meff)*(max(out_param.abstol, out_param.reltol*abs(q + errest(meff)))...
        - max(out_param.abstol, out_param.reltol*abs(q - errest(meff))))/...
        (max(out_param.abstol, out_param.reltol*abs(q + errest(meff)))...
        + max(out_param.abstol, out_param.reltol*abs(q - errest(meff)))); % Optimal estimator
   appxinteg(meff)=q;

   if 4*errest(meff)^2/(max(out_param.abstol, out_param.reltol*abs(q + errest(meff)))...
        + max(out_param.abstol, out_param.reltol*abs(q - errest(meff))))^2 <= 1
      out_param.time=toc(t_start);
      is_done = true;
   elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
      out_param.exit(1) = true;
   end
end

% Decode the exit structure
exit_str=2.^(0:exit_len-1).*out_param.exit;
exit_str(out_param.exit==0)=[];
if numel(exit_str)==0;
    out_param.exitflag=0;
else
    out_param.exitflag=exit_str;
end

out_param = rmfield(out_param,'exit');

out_param.time=toc(t_start);
end


%% Parsing for the input of cubLattice_g
function [f,hyperbox, out_param] = cubLattice_g_param(r_lag,varargin)

% Default parameter values
default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
default.measure  = 'uniform';
default.transf = 1;% default transformation (box-to-ball or box-to-sphere)
default.radius = 1;% radius of the ball or sphere
default.abstol  = 1e-4;
default.reltol  = 1e-2;
default.shift  = rand;
default.mmin  = 10;
default.mmax  = 20;
default.fudge = @(m) 5*2.^-m;
default.transform = 'Baker';

if numel(varargin)<2
    help cubLattice_g
    warning('GAIL:cubLattice_g:fdnotgiven',...
        'At least, function f and hyperbox need to be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    hyperbox = default.hyperbox;
else
    f = varargin{1};
    if ~gail.isfcn(f)
        warning('GAIL:cubLattice_g:fnotfcn',...
            'The given input f was not a function. Example for f(x)=x^2:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
    elseif numel(varargin) == 2
        out_param.f=f;
        hyperbox = varargin{2};
        if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(size(hyperbox,2)<601)
            warning('GAIL:cubLattice_g:hyperbox_error1',...
                'The hyperbox must be a real matrix of size 2xd where d can not be greater than 600. Example for f(x)=x^2:')
            f = @(x) x.^2;
            out_param.f=f;
            hyperbox = default.hyperbox;
        end
     else
        out_param.f = f;
        hyperbox = varargin{2}; % hyperbox validation will be done above
        default.shift = rand(1, size(hyperbox, 2)); % The shift needs to take the dimension of the problem
    end
end

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin(3:end);
    for j=1:numel(varargin)-2
    validvarargin=validvarargin && (isnumeric(in3{j}) ...
        || ischar(in3{j}) || isstruct(in3{j}) || gail.isfcn(in3{j}));
    end
    if ~validvarargin
        warning('GAIL:cubLattice_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in3=varargin{3};
end

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
  f_addParamVal = @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin
    out_param.measure = default.measure;
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.shift = default.shift;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;
    out_param.fudge = default.fudge;
    out_param.transform = default.transform;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'hyperbox',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal','uniform ball',...
            'uniform ball_box','uniform ball_normal','uniform sphere',...
            'uniform sphere_box','uniform sphere_normal'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'shift',default.shift,@isnumeric);
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
        addOptional(p,'transform',default.transform,...
            @(x) any(validatestring(x, {'id','Baker','C0','C1','C1sin'})));
    else
        if isstruct(in3) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal','uniform ball',...
            'uniform ball_box','uniform ball_normal','uniform sphere',...
            'uniform sphere_box','uniform sphere_normal'})));
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'shift',default.shift,@isnumeric);
        f_addParamVal(p,'mmin',default.mmin,@isnumeric);
        f_addParamVal(p,'mmax',default.mmax,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@gail.isfcn);
        f_addParamVal(p,'transform',default.transform,...
            @(x) any(validatestring(x, {'id','Baker','C0','C1','C1sin'})));
    end
    parse(p,f,hyperbox,varargin{3:end});
    out_param = p.Results;
end

% Force measure to be one of the allowed ones
if ~(strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') || ...
        strcmp(out_param.measure,'uniform ball') || ...
        strcmp(out_param.measure,'uniform ball_box') || ...
        strcmp(out_param.measure,'uniform ball_normal') || ...
        strcmp(out_param.measure,'uniform sphere') || ...
        strcmp(out_param.measure,'uniform sphere_box') || ...
        strcmp(out_param.measure,'uniform sphere_normal'))
    warning('GAIL:cubLattice_g:notmeasure',['Given measure is not allowed.' ...
            ' Using default measure ' num2str(default.measure)])
    out_param.measure = default.measure;
end

% simplifying out_param.measure and storing which transformation should be
% applied
if strcmp(out_param.measure,'uniform ball') || strcmp(out_param.measure,'uniform sphere')
    out_param.transf = default.transf;
elseif strcmp(out_param.measure,'uniform ball_box')
    out_param.transf = 1; % one means from a box
    out_param.measure = 'uniform ball';
elseif strcmp(out_param.measure,'uniform ball_normal')
    out_param.transf = 2; % two means from normal
    out_param.measure = 'uniform ball';
elseif strcmp(out_param.measure,'uniform sphere_box')
    out_param.transf = 1;
    out_param.measure = 'uniform sphere';
elseif strcmp(out_param.measure,'uniform sphere_normal')
    out_param.transf = 2;
    out_param.measure = 'uniform sphere';
end

%validating hyperbox
if strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') % 'uniform box' or 'normal box'
    out_param.d = size(hyperbox,2);
    if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(out_param.d<1111)
        warning('GAIL:cubLattice_g:hyperbox_error2',...
            'When measure is ''uniform'' or ''normal'', the hyperbox must be a real matrix of size 2 x d where d can not be greater than 1111. Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end
else % 'uniform ball' or 'uniform sphere'
    out_param.d = size(hyperbox,2);
    if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==1) || ~(out_param.d<1112) % size(hyperbox,2) is actually equal to d+1 (the extra value is the radius)
        warning('GAIL:cubLattice_g:hyperbox_error3',...
            'When measure is ''uniform ball'' or ''uniform sphere'', the hyperbox must be a real matrix of size 1 x (d+1) where d can not be greater than 1111.  Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end

    out_param.radius = hyperbox(out_param.d);
    hyperbox = hyperbox(:,1:out_param.d-1); % removing the last value is the radius, which is the radius
    out_param.d = out_param.d - 1; % storing the rigth dimension of the ball or sphere
    out_param.shift = out_param.shift(1:end-1);

    if strcmp(out_param.measure,'uniform ball') && out_param.d <= 0
        warning('GAIL:cubLattice_g:dimensionequalszero',...
            'When measure is ''uniform ball'', the number of dimentions must be a positive values.  Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end

    if strcmp(out_param.measure,'uniform sphere') && out_param.d <= 1
        warning('GAIL:cubLattice_g:dimensionsmallerthan2',...
            'When measure is ''uniform sphere'', the number of dimentions must be at least 2.  Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end

    if ~isfinite(out_param.radius) || out_param.radius <= 0.0
        warning('GAIL:cubLattice_g:infiniteradius',...
            'When measure is ''uniform ball'' or ''uniform sphere'', the radius must a finite positive real number. Default value for the radius will be used:')
        out_param.radius = default.radius;
    end

    if strcmp(out_param.measure,'uniform sphere') && out_param.transf == 1 % box-to-sphere transformation
        % setting out_param.d to be the dimension of the box over which the
        % integral will actually be computed
        out_param.d = out_param.d - 1;
        out_param.shift = out_param.shift(1:end-1);
    end
end

% Force absolute tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:cubLattice_g:abstolnonpos',['Absolute tolerance cannot be negative.' ...
            ' Using default absolute tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol < 0) || (out_param.reltol > 1)
    warning('GAIL:cubLattice_g:reltolnonunit',['Relative tolerance should be chosen in [0,1].' ...
            ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Checks if shift is a vector in [0,1)^d
if ~(all(out_param.shift < 1) && all(out_param.shift >= 0) && size(out_param.shift, 2) == out_param.d)
    warning('GAIL:cubLattice_g:shift',['The shift should be a vector ' ...
            'of size 1 x d in [0,1)^d.' ...
            ' Using default shift ' num2str(default.shift)])
    out_param.shift = default.shift;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('GAIL:cubLattice_g:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force mmin to be integer greater than r_lag (so that l_star=mmin-r_lag>=0)
if out_param.mmin < r_lag
    warning('GAIL:cubLattice_g:lowmminrlag',['The minimum starting exponent ' ...
            'should be at least ' num2str(r_lag) '.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 20
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin)
    warning('GAIL:cubLattice_g:wrongmmax',['The maximum exponent for the budget should be an integer bigger than mmin and smaller than the allowed by gail.lattice_gen.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('GAIL:cubLattice_g:fudgenonpos',['The fudge factor should be a positive function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force transform to only be id, Baker, C0, C1 or C1sin
if ~(strcmp(out_param.transform,'id') || strcmp(out_param.transform,'Baker') || strcmp(out_param.transform,'C0') || strcmp(out_param.transform,'C1') || strcmp(out_param.transform,'C1sin') )
    warning('GAIL:cubLattice_g:notmeasure',['The periodizing transformations can only be id, Baker, C0, C1 or C1sin.' ...
            ' Using default error tolerance ' num2str(default.transform)])
    out_param.transform = default.transform;
end

% Checking on pure absolute/relative error
if (out_param.abstol==0) && (out_param.reltol==0)
    warning('GAIL:cubLattice_g:tolzeros',['Absolute and relative error tolerances can not be simultaniusly 0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ' and relative tolerance ' num2str(default.reltol)])
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
end

% Checking on the hyperbox given the measure
if (strcmp(out_param.measure,'uniform')) && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubLattice_g:hyperboxnotfinite',['If uniform measure, hyperbox must be of finite volume.' ...
            ' Using default hyperbox:'])
    disp([zeros(1,out_param.d);ones(1,out_param.d)])
    hyperbox = [zeros(1,out_param.d);ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(any(isfinite(hyperbox)))>0)
    warning('GAIL:cubLattice_g:hyperboxfinite',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(hyperbox(1,:)==hyperbox(2,:)) || any(hyperbox(1,:)>hyperbox(2,:)))
    warning('GAIL:cubLattice_g:hyperboxnormalwrong',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'uniform ball') || strcmp(out_param.measure,'uniform sphere'))...
        && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubLattice_g:infinitecoordinateforthecenter',['If uniform ball or sphere measure, all the coordinates of the center must be finite.' ...
            ' Using the origin as the center:'])
    % out_param.d should not be used here because this variable stores the
    % dimension of the box over which the integral will actually be
    % computed, whih may be different from the dimesion of the sphere
    hyperbox = zeros(1,size(hyperbox,2));
end

end
