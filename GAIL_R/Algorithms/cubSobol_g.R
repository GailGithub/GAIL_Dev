library(randtoolbox)
library(pracma)
cubSobol_g = function(hyperbox = c(0,1), measure = 'uniform', 
                      abstol = 1e-4, reltol = 1e-2, mmin = 10, mmax = 24, fudge = function(m) {5*2^(-m)}, toltype = 'max', theta = 1){
#function [q,out_param] = cubSobol_g(varargin)
#cubSobol_G Quasi-Monte Carlo method using Sobol' cubature over the
#d-dimensional region to integrate within a specified generalized error
#tolerance with guarantees under Walsh-Fourier coefficients cone decay
#assumptions
#
#   [q,out_param] = CUBSOBOL_G(f,hyperbox) estimates the integral of f
#   over the d-dimensional region described by hyperbox, and with an error
#   guaranteed not to be greater than a specific generalized error tolerance,
#   tolfun:=max(abstol,reltol*| integral(f) |). Input f is a function handle. f should
#   accept an n x d matrix input, where d is the dimension and n is the 
#   number of points being evaluated simultaneously. The input hyperbox is
#   a 2 x d matrix, where the first row corresponds to the lower limits 
#   and the second row corresponds to the upper limits of the integral.
#   Given the construction of Sobol' sequences, d must be a positive 
#   integer with 1<=d<=1111.
#
#   q = CUBSOBOL_G(f,hyperbox,measure,abstol,reltol)
#   estimates the integral of f over the hyperbox. The answer
#   is given within the generalized error tolerance tolfun. All parameters
#   should be input in the order specified above. If an input is not specified,
#   the default value is used. Note that if an input is not specified,
#   the remaining tail cannot be specified either. Inputs f and hyperbox 
#   are required. The other optional inputs are in the correct order:
#   measure,abstol,reltol,mmin,mmax,fudge,toltype and
#   theta.
#
#   q = CUBSOBOL_G(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
#   estimates the integral of f over the hyperbox. The answer
#   is given within the generalized error tolerance tolfun. All the field-value
#   pairs are optional and can be supplied in any order. If an input is not
#   specified, the default value is used.
#
#   q = CUBSOBOL_G(f,hyperbox,in_param) estimates the integral of f over the
#   hyperbox. The answer is given within the generalized error tolerance tolfun.
# 
#   Input Arguments
#
#     f --- the integrand whose input should be a matrix n x d where n is
#     the number of data points and d the dimension, which cannot be
#     greater than 1111. By default f is f=@ x.^2.
#
#     hyperbox --- the integration region defined by its bounds. It must be
#     a 2 x d matrix, where the first row corresponds to the lower limits 
#     and the second row corresponds to the upper limits of the integral.
#     The default value is [0;1].
#
#     in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
#     measure of a uniformly distributed random variable in the hyperbox
#     or normally distributed with covariance matrix I_d. The only possible
#     values are 'uniform' or 'normal'. For 'uniform', the hyperbox must be
#     a finite volume while for 'normal', the hyperbox can only be defined as 
#     (-Inf,Inf)^d. By default it is 'uniform'.
#
#     in_param.abstol --- the absolute error tolerance, abstol>=0. By 
#     default it is 1e-4.
#
#     in_param.reltol --- the relative error tolerance, which should be
#     in [0,1]. Default value is 1e-2.
# 
#   Optional Input Arguments
# 
#     in_param.mmin --- the minimum number of points to start is 2^mmin.
#     The cone condition on the Fourier coefficients decay requires a
#     minimum number of points to start. The advice is to consider at least
#     mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
#     default it is 10.
#
#     in_param.mmax --- the maximum budget is 2^mmax. By construction of
#     the Sobol' generator, mmax is a positive integer such that
#     mmin<=mmax<=53. The default value is 24.
#
#     in_param.fudge --- the positive function multiplying the finite 
#     sum of Fast Walsh Fourier coefficients specified in the cone of functions.
#     This input is a function handle. The fudge should accept an array of
#     nonnegative integers being evaluated simultaneously. For more
#     technical information about this parameter, refer to the references.
#     By default it is @(m) 5*2.^-m.
#
#     in_param.toltype --- this is the generalized tolerance function.
#     There are two choices, 'max' which takes
#     max(abstol,reltol*| integral(f) | ) and 'comb' which is the linear combination
#     theta*abstol+(1-theta)*reltol*| integral(f) | . Theta is another 
#     parameter to be specified with 'comb'(see below). For pure absolute
#     error, either choose 'max' and set reltol = 0 or choose 'comb' and set
#     theta = 1. For pure relative error, either choose 'max' and set 
#     abstol = 0 or choose 'comb' and set theta = 0. Note that with 'max',
#     the user can not input abstol = reltol = 0 and with 'comb', if theta = 1
#     abstol con not be 0 while if theta = 0, reltol can not be 0.
#     By default toltype is 'max'.
# 
#     in_param.theta --- this input is parametrizing the toltype 
#     'comb'. Thus, it is only active when the toltype
#     chosen is 'comb'. It establishes the linear combination weight
#     between the absolute and relative tolerances
#     theta*abstol+(1-theta)*reltol*| integral(f) |. Note that for theta = 1, 
#     we have pure absolute tolerance while for theta = 0, we have pure 
#     relative tolerance. By default, theta=1.
#
#   Output Arguments
#
#     q --- the estimated value of the integral.
#
#     out_param.d --- dimension over which the algorithm integrated.
#
#     out_param.n --- number of Sobol' points used for computing the
#     integral of f.
#
#     out_param.bound_err --- predicted bound on the error based on the cone
#     condition. If the function lies in the cone, the real error will be
#     smaller than generalized tolerance.
#
#     out_param.time --- time elapsed in seconds when calling cubSobol_g.
#
#     out_param.exitflag --- this is a binary vector stating whether
#     warning flags arise. These flags tell about which conditions make the
#     final result certainly not guaranteed. One flag is considered arisen
#     when its value is 1. The following list explains the flags in the
#     respective vector order:
#
#                      1 : If reaching overbudget. It states whether
#                       the max budget is attained without reaching the
#                       guaranteed error tolerance.
#      
#                       2 : If the function lies outside the cone. In
#                       this case, results are not guaranteed. For more
#                       information about the cone definition, check the
#                       article mentioned below.
# 
#  Guarantee
# This algorithm computes the integral of real valued functions in [0,1)^d
# with a prescribed generalized error tolerance. The Walsh-Fourier
# coefficients of the integrand are assumed to be absolutely convergent. If
# the algorithm terminates without warning messages, the output is given
# with guarantees under the assumption that the integrand lies inside a
# cone of functions. The guarantee is based on the decay rate of the
# Walsh-Fourier coefficients. For more details on how the cone is defined,
# please refer to the references below.
# 
#  Examples
# 
# Example 1:
# % Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:
# % 
# >> f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)]; 
# >> q = cubSobol_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
# >> check = abs(exactsol-q) < 1e-5
# check = 1
# 
# 
# Example 2:
# % Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2
# in the interval R^3 where x1, x2 and x3 are normally distributed:
# % 
# >> f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
# >> q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3); exactsol = 1;
# >> check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max')
# check = 1
# 
# 
# Example 3: 
# % Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
# interval [-1,2)^2:
# % 
# >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
# >> q = cubSobol_g(f,hyperbox,'uniform',1e-3,1e-2); exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
# >> check = abs(exactsol-q) < gail.tolfun(1e-3,1e-2,1,exactsol,'max')
# check = 1
#
#
# Example 4: 
# % Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
# sigma=0.05 and T=1.
# 
# >> f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [-inf(1,1);inf(1,1)];
# >> q = cubSobol_g(f,hyperbox,'normal',1e-4,1e-2); price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
# >> check = abs(price-q) < gail.tolfun(1e-4,1e-2,1,price,'max')
# check = 1
#
#
# Example 5:
# % Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the interval
# [0,1)^5 with pure absolute error 1e-5.
# 
# >> f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
# >> q = cubSobol_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
# >> check = abs(exactsol-q) < 1e-5
# check = 1
#
#
#   See also CUBLATTICE_G, CUBMC_G, MEANMC_G, MEANMCBER_G, INTEGRAL_G
# 
#  References
#
#   [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable adaptive
#   cubature using digital sequences," 2014. Submitted for publication:
# %   arXiv:1410.8615.
#
#   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
#   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
#   GAIL: Guaranteed Automatic Integration Library (Version 2.1)
#   [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
# %
#   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
#   Research via Supportable Scientific Software," Journal of Open Research
#   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
#
#   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
#   Mathematical Software" [Course Slides], Illinois Institute of
#   Technology, Chicago, IL, 2013. Available from
#   http://code.google.com/p/gail/ 
# %
#   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
#   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
#   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
#   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
#   Workshop On Sustainable Software for Science: Practice And Experiences
#   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
#   pp. 1-21, 2014.
#
#   If you find GAIL helpful in your work, please support us by citing the
#   above papers, software, and materials.
#
#   Authors: Anthony Karahalios and Luana Terra
#
##Initial important cone factors and Check-initialize parameters
r_lag = 4; #distance between coefficients summed and those computed
cubSobol_g_param() #Get the parameters
l_star = my_list$"out_param"["mmin"] - r_lag; # Minimum gathering of points for the sums of DFWT
  
if (out_param.measure == 'normal') {
f = function(x) {qnorm(x)^2}
}
else if (out_param.measure == 'uniform') {
Cnorm = prod(hyperbox[,2]-hyperbox[,1]      # # # # # hyperbox needs to be two columns instead of two rows in R
f = function(x) {Cnorm *(hyperbox[,1]+(hyperbox[,2]-hyperbox[,1])*x)^2}
}

##Main algorithm
Stilde=rep(0,out_param.mmax-out_param.mmin+1); #initialize sum of DFWT terms
Cstilde_low=matrix(-Inf,1,out_param.mmax-l_star+1); #initialize #various sums of DFWT terms for necessary conditions
CStilde_up=matrix(Inf,1,out_param.mmax-l_star+1); #initialize various #sums of DFWT terms for necessary conditions
errest=rep(0,out_param.mmax-out_param.mmin+1); #initialize error estimates
appxinteg=rep(0,out_param.mmax-out_param.mmin+1); #initialize approximations to integral
exit_len = 2;
#out_param.exit=false(1,exit_len); #we start the algorithm with #all warning flags down

##Initial points and FWT
out_param.n=2^out_param.mmin; #total number of points to start with
n0=out_param.n; #initial number of points
xpts=sobol(n0,out_param.d,scrambling=1)
####xpts=sobstr(1:n0,1:out_param.d); %grab Sobol' points     # # # # # # We are trying to use xpts directly without sobstr
y=f(xpts); #evaluate integrand
yval=y;
}

## Compute initial FWT
for (l in 0:out_param.mmin-1){
  nl=2^l;
  nmminlm1=2^(out_param.mmin-l-1);
  ptind=matrix(rep(matrix(c(rep(TRUE,nl),rep(FALSE,nl))),nmminlm1));
  evenval=y[ptind];
  oddval=y[!ptind];
  y[ptind]=(evenval+oddval)/2;
  y[!ptind]=(evenval-oddval)/2;
}
#y now contains the FWT coefficients

## Approximate integral
q = mean(yval);
appxinteg[2] = q;

## Create kappanumap implicitly from the data
kappanumap = (1:out_param.n); #initialize map
for (l in seq((out_param.mmin-1),1,-1)) {
  nl = 2^l;
  oldone = abs(y[kappanumap[(2:nl)]]); #earlier values of kappa, don't touch first ones
  newone = abs(y[kappanumap[(nl+2):(2*nl)]]); #later values of kappa,
  flip = which(newone>oldone); #which in the pair are the larger ones
  if (length(flip)!=0) {
    flipall = bsxfun("+",flip,seq(0,2^(out_param.mmin-1),2^(l+1)));
    flipall=as.vector(flipall);
    temp=kappanumap[nl+1+flipall]; #then flip
    kappanumap[nl+1+flipall]=kappanumap[1+flipall]; #them
    kappanumap[1+flipall]=temp; #around
  }
}
## Compute Stilde
nllstart=2^(out_param.mmin-r_lag-1);  #no good equiv for int64 w/o package
Stilde[1]=sum(abs(y[kappanumap[nllstart+1:2*nllstart]]));
out_param.bound_err=out_param.fudge(out_param.mmin)*Stilde[1];
errest[1]=out_param.bound_err;

# Necessary conditions
for (l in l_star:out_param.mmin) { # Storing the information for the necessary conditions
C_low = (1+out_param.fudge(out_param.mmin-l))/(1+2*out_param.fudge(out_param.mmin-l));
C_up = (1+out_param.fudge(out_param.mmin-l));
CStilde_low[l-l_star+1] = C_low*sum(abs(y[kappanumap[2^(l-1)+1:2^l]]));
CStilde_up[l-l_star+1] = C_up*sum(abs(y[kappanumap[2^(l-1)+1:2^l]]));
}
###if (any(CStilde_low(:) > CStilde_up(:))) { ####Unsure about this part
###out_param.exit(2) = true; #####Unsure about this part
}

##Defining the Parameters
cubSobol_g_param = function(r_lag, hyperbox = c(0,1), measure = 'uniform', 
                            abstol = 1e-4, reltol = 1e-2, mmin = 10, mmax = 24, fudge = function(m) {5*2^(-m)}, toltype = 'max', theta = 1){ 
  #To define the variables but we still need to code for the checks****
  #  default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
  #  default.measure  = 'uniform';
  #  default.abstol  = 1e-4;
  #  default.reltol  = 1e-2;
  #  default.mmin  = 10;
  #  default.mmax  = 24;
  #  default.fudge = @(m) 5*2.^-m;
  #  default.toltype  = 'max';
  #  default.theta  = 1;
  out_param.measure = measure;
  out_param.abstol = abstol;
  out_param.reltol = reltol;
  out_param.mmin = mmin;
  out_param.mmax = mmax;  
  out_param.fudge = fudge;
  out_param.toltype = toltype;
  out_param.theta = theta;
  out_param=c("measure"=out_param.measure,"abstol"=out_param.abstol,"reltol"=out_param.reltol,"mmin"=out_param.mmin,"mmax"=out_param.mmax,
              "fudge"=out_param.fudge,"toltype"=out_param.toltype,"theta"=out_param.theta)
  out_list = list("out_param"=out_param,"f"=f,"hyperbox"=hyperbox)
  return(out_list)
}