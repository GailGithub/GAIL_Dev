#MEANMC_G Monte Carlo method to estimate the mean of a random variable
#
#   tmu = MEANMC_G(Yrand) estimates the mean, mu, of a random variable Y to
#   within a specified generalized error tolerance, 
#   tolfun:=max(abstol,reltol*| mu |), i.e., | mu - tmu | <= tolfun with
#   probability at least 1-alpha, where abstol is the absolute error
#   tolerance, and reltol is the relative error tolerance. Usually the
#   reltol determines the accuracy of the estimation, however, if the | mu |
#   is rather small, the abstol determines the accuracy of the estimation.
#   The default values are abstol=1e-2, reltol=1e-1, and alpha=1%. Input
#   Yrand is a function handle that accepts a positive integer input n and
#   returns an n x 1 vector of IID instances of the random variable Y.
#
#   tmu = MEANMC_G(Yrand,abstol,reltol,alpha) estimates the mean of a
#   random variable Y to within a specified generalized error tolerance
#   tolfun with guaranteed confidence level 1-alpha using all ordered
#   parsing inputs abstol, reltol, alpha.
#   
#   tmu = MEANMC_G(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha)
#   estimates the mean of a random variable Y to within a specified
#   generalized error tolerance tolfun with guaranteed confidence level
#   1-alpha. All the field-value pairs are optional and can be supplied in
#   different order, if a field is not supplied, the default value is used.
#
#   [tmu, out_param] = MEANMC_G(Yrand,in_param) estimates the mean of a
#   random variable Y to within a specified generalized error tolerance
#   tolfun with the given parameters in_param and produce the estimated
#   mean tmu and output parameters out_param. If a field is not specified,
#   the default value is used.
#
#   Input Arguments
#
#     Yrand --- the function for generating n IID instances of a random
#     variable Y whose mean we want to estimate. Y is often defined as a
#     function of some random variable X with a simple distribution. The
#     input of Yrand should be the number of random variables n, the output
#     of Yrand should be n function values. For example, if Y = X.^2 where X
#     is a standard uniform random variable, then one may define Yrand =
#     @(n) rand(n,1).^2.
#
#     in_param.abstol --- the absolute error tolerance, which should be
#     positive, default value is 1e-2.
#
#     in_param.reltol --- the relative error tolerance, which should be
#    between 0 and 1, default value is 1e-1.
#
#     in_param.alpha --- the uncertainty, which should be a small positive
#     percentage. default value is 1%.
#
#   Optional Input Arguments
#
#     in_param.fudge --- standard deviation inflation factor, which should
#     be larger than 1, default value is 1.2.
#
#     in_param.nSig --- initial sample size for estimating the sample
#     variance, which should be a moderate large integer at least 30, the
#     default value is 1e4.
#
#     in_param.n1 --- initial sample size for estimating the sample mean,
#     which should be a moderate large positive integer at least 30, the
#     default value is 1e4.
#
#     in_param.tbudget --- the time budget in seconds to do the two-stage
#     estimation, which should be positive, the default value is 100 seconds.
#
#     in_param.nbudget --- the sample budget to do the two-stage
#     estimation, which should be a large positive integer, the default
#     value is 1e9.
#
#   Output Arguments
#
#     tmu --- the estimated mean of Y.
#
#     out_param.tau --- the iteration step.
#
#     out_param.n --- the sample size used in each iteration.
#
#     out_param.nremain --- the remaining sample budget to estimate mu. It was
#     calculated by the sample left and time left.
#
#     out_param.ntot --- total sample used.
#
#     out_param.hmu --- estimated mean in each iteration.
#
#     out_param.tol --- the reliable upper bound on error for each iteration.
#
#     out_param.var --- the sample variance.
#
#     out_param.exit --- the state of program when exiting.
#       
#                      0   Success
#      
#                      1   Not enough samples to estimate the mean
#
#     out_param.kurtmax --- the upper bound on modified kurtosis.
#
#     out_param.time --- the time elapsed in seconds.
#
#     out_param.flag --- parameter checking status
#      
#                           1  checked by meanMC_g
#
#  Guarantee
# This algorithm attempts to calculate the mean, mu, of a random variable
# to a prescribed error tolerance, tolfun:= max(abstol,reltol*| mu |), with
# guaranteed confidence level 1-alpha. If the algorithm terminated without
# showing any warning messages and provide an answer tmu, then the follow
# inequality would be satisfied:
# 
# Pr(| mu - tmu | <= tolfun) >= 1-alpha
#
# The cost of the algorithm, N_tot, is also bounded above by N_up, which is
# defined in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And
# the following inequality holds:
# 
# Pr (N_tot <= N_up) >= 1-beta
#
# Please refer to our paper for detailed arguments and proofs.
#
# Examples
#
# Example 1:
# If no parameters are parsed, help text will show up as follows:
# >> meanMC_g
# ***Monte Carlo method to estimate***
#
#
# Example 2:
# Calculate the mean of x^2 when x is uniformly distributed in
# [0 1], with the absolute error tolerance = 1e-3 and uncertainty 5%.
#
# >> in_param.reltol=0; in_param.abstol = 1e-3;
# >> in_param.alpha = 0.05; Yrand=@(n) rand(n,1).^2;
# >> tmu=meanMC_g(Yrand,in_param)
# tmu = 0.33***
#
#
# Example 3:
# Calculate the mean of exp(x) when x is uniformly distributed in
# [0 1], with the absolute error tolerance 1e-3.
#
# >> tmu=meanMC_g(@(n)exp(rand(n,1)),1e-3,0)
# tmu = 1.71***
#
#
# Example 4:
# Calculate the mean of cos(x) when x is uniformly distributed in
# [0 1], with the relative error tolerance 1e-2 and uncertainty 0.05.
#
# >> tmu=meanMC_g(@(n)cos(rand(n,1)),'reltol',1e-2,'abstol',0,'alpha',0.05)
# tmu = 0.84***
#
#
#   See also FUNAPPX_G, INTEGRAL_G, CUBMC_G, MEANMCBER_G, CUBSOBOL_G, CUBLATTICE_G
#
#  References
#
#   [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, "Guaranteed
#   conservative fixed width confidence intervals via Monte Carlo
#   sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
#   Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin,
#   2014. arXiv:1208.4318 [math.ST]
#            
#   [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
#   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
#   GAIL: Guaranteed Automatic Integration Library (Version 2.1) [MATLAB
#   Software], 2015. Available from http://gailgithub.github.io/GAIL_Dev/
#
#   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
#   Research via Supportable Scientific Software," Journal of Open Research
#   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
#
#   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
#   Mathematical Software" [Course Slides], Illinois Institute of
#   Technology, Chicago, IL, 2013. Available from
#   http://gailgithub.github.io/GAIL_Dev/ 
#
#   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
#   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
#   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
#   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
#   Workshop On Sustainable Software for Science: Practice and Experiences
#   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
#   pp. 1-21, 2014.
#
#   If you find GAIL helpful in your work, please support us by citing the
#   above papers, software, and materials.
#
#Create function with default values given
meanMC_g = function(Yrand = function(n) {runif(n)^2},abstol = 1e-2, reltol = 1e-1, nSig = 1e4, n1 = 1e4, alpha = 0.01, fudge = 1.2, tbudget = 100, nbudget = 1e9)
{
  tstart = proc.time(); #start the clock
  ####[Yrand, out_param] = meanMC_g_param(varargin{:})
  meanMC_g_param();
  ####out_param = c(abstol,reltol,nSig,n1,alpha,fudge,tbucket,nbucket)
  
  n1 = 2 #let it run once to load all the data. warm up the machine.
  Yrand(n1);
  nsofar = n1;
  
  ntry = 10; #try several samples to get the time
  tstart1 = proc.time()
  Yrand(ntry); 
  ttry = proc.time() - tstart1;
  tpern = ttry[3]/ntry; # calculate time per sample
  nsofar = nsofar+ntry; # update n so far
  out_param.exit = 0;
  
  if(tpern<1e-7); #each sample use very very little time
  booster = 8;
  tstart2 = proc.time();Yrand(ntry*booster);ttry2 = proc.time() - tstart2;
  ntry = ntry*c(1,booster);
  ttry = c(ttry,ttry2);# take eight times more samples to try
  else if(tpern>=1e-7 && tpern<1e-5) #each sample use very little time
    booster = 6;
  tstart3 = proc.time();Yrand(ntry*booster);ttry2 = proc.time() - tstart3;
  ntry = ntry*c(1,booster);
  ttry = c(ttry,ttry2);# take six times more samples to try    
  else if(tpern>=1e-5 && tpern<1e-3) #each sample use little time
    booster = 4;
  tstart4 = proc.time();Yrand(ntry*booster);ttry2 = proc.time() - tstart4;
  ntry = ntry*c(1,booster);
  ttry = c(ttry,ttry2);# take four times more samples to try
  else if(tpern>=1e-3 && tpern<1e-1) #each sample use moderate time
    booster = 2;
  tstart5 = proc.time();Yrand(ntry*booster);ttry2 = proc.time() - tstart5;
  ntry = ntry*c(1,booster);
  ttry = c(ttry,ttry2);# take two times more samples to try
  else #each sample use lots of time, stop try
    end()
  #####[tmu,out_param] =  meanmctolfun(Yrand,out_param,ntry,ttry,nsofar,tstart);
  #####end
  
meanmctolfun = function(Yrand(),out_param,ntry,ttry,nsofar,tstart) {
  tstart6 = proc.time();
  Yval = Yrand(out_param[4]);# get samples to estimate variance 
  t_sig = proc.time()-tstart6;#get the time for calculating nSig function values.
  nsofar = nsofar+out_param[4];# update the samples that have been used
  source('estsamplebudget.R') #File needs to be in Working Directory
  out_param.nremain = estsamplebudget(out_param[8],
  out_param[9],c(ntry,out_param[4]),nsofar,tstart,c(ttry[3],t_sig[3]);
#update the nremain could afford until now
out_param.var = var(Yval);# calculate the sample variance--stage 1
sig0 = sqrt(out_param.var);# standard deviation
sig0up = out_param[7].*sig0;# upper bound on the standard deviation
alpha_sig = out_param[6]/2;# the uncertainty for variance estimation
alphai = (out_param[6]-alpha_sig)/(2*(1-alpha_sig));
#uncertainty to do iteration
out_param.kurtmax = (out_param[4]-3)/(out_param[4]-1) ...
+ ((alpha_sig*out_param[4])/(1-alpha_sig))...
*(1-1/out_param[7]^2)^2;
#the upper bound on the modified kurtosis
npcmax = 1e6;#constant to do iteration and mean calculation4
if(out_param[3]==0) {
alphai = 1-(1-out_param[6])/(1-alpha_sig)
  if(sig0up == 0) { # if the variance is zero, just take n_sigma samples
out_param.n = out_param[4]
  }
  else {
  toloversig = out_param.abstol/sig0up;
# absolute error tolerance over sigma
out_param.n = nchebe(toloversig,alphai,out_param.kurtmax);
if(out_param.n > out_param.nremain) {
out_param.exit=1; #pass a flag
meanMC_g_err(out_param); # print warning message
out_param.n = out_param.nremain;# update n
}
}
source('evalmean.R') #must be in Working Directory
tmu = evalmean(Yrand,out_param.n,npcmax);#evaluate the mean
nsofar = nsofar+out_param.n;#total sample used
else {
  alphai = (out_param[6]-alpha_sig)/(2*(1-alpha_sig));
#uncertainty to do iteration
eps1 = ncbinv(out_param[5],alphai,out_param.kurtmax);
#tolerance for initial estimation
out_param.tol = c(sig0up*eps1);
#the width of initial confidence interval for the mean
i=1;
out_param.n[i] = out_param[5];# initial sample size to do iteration
while(true) {
out_param.tau = i;#step of the iteration
if(out_param.n[i] > out_param.nremain) {
# if the sample size used for initial estimation is
# larger than nremain, print warning message and use nremain
out_param.exit=1; #pass a flag
meanMC_g_err(out_param); # print warning message
out_param.n[i] = out_param.nremain;# update n
tmu = evalmean(Yrand,out_param.n(i),npcmax);#evaluate the mean
nsofar = nsofar+out_param.n[i];#total sample used
break;
}
out_param.hmu = c()
out_param.hmu[i] = evalmean(Yrand,out_param.n[i],npcmax);#evaluate mean
nsofar = nsofar+out_param.n[i];
out_param.nremain = out_param.nremain-out_param.n[i];#update n so far and nremain
errtype = 'max';
# error type, see the function 'tolfun' at Algoithms/+gail/ directory
# for more info
theta  = 0;# relative error case
deltaplus = (gail.tolfun(out_param.abstol,out_param.reltol,...
                         theta,out_param.hmu(i) - out_param.tol(i),errtype)...
             +gail.tolfun(out_param.abstol,out_param.reltol,...
                          theta,out_param.hmu(i) + out_param.tol(i),errtype))/2;
# a combination of tolfun, which used to decide stopping time
if deltaplus >= out_param.tol(i) # stopping criterion
deltaminus= (gail.tolfun(out_param.abstol,out_param.reltol,...
                         theta,out_param.hmu(i) - out_param.tol(i),errtype)...
             -gail.tolfun(out_param.abstol,out_param.reltol,...
                          theta,out_param.hmu(i) + out_param.tol(i),errtype))/2;
# the other combination of tolfun, which adjust the hmu a bit
tmu = out_param.hmu(i)+deltaminus;
break;
else
  out_param.tol(i+1) = min(out_param.tol(i)/2,max(out_param.abstol,...
                                                  0.95*out_param.reltol*abs(out_param.hmu(i))));
i=i+1;
end
toloversig = out_param.tol(i)/sig0up;#next tolerance over sigma
alphai = (out_param.alpha-alpha_sig)/(1-alpha_sig)*2.^(-i);
#update the next uncertainty
out_param.n(i) = nchebe(toloversig,alphai,out_param.kurtmax);
end
end
#get the next sample size needed
out_param.ntot = nsofar;%total sample size used
out_param.time=toc(tstart); %elapsed time
}

ncbinv = function(n1,alpha1,output_param.kurtmax) {
#This function calculate the reliable upper bound on error when given
#Chebyshev and Berry-Esseen inequality and sample size n.
NCheb_inv = 1/sqrt(n1*alpha1);
#use Chebyshev inequality
A=18.1139;
A1=0.3328;
A2=0.429; # three constants in Berry-Esseen inequality
M3upper=output_param.kurtmax^(3/4);
#using Jensen's inequality to bound the third moment
BEfun=@(logsqrtb)gail.stdnormcdf(n1.*logsqrtb)...
    +min(A1*(M3upper+A2), ...
    A*M3upper./(1+(sqrt(n1).*logsqrtb).^3))/sqrt(n1)...
    - alpha1/2;
# Berry-Esseen inequality
logsqrtb_clt=log(sqrt(gail.stdnorminv(1-alpha1/2)/sqrt(n1)));
#use CLT to get tolerance
NBE_inv = exp(2*fzero(BEfun,logsqrtb_clt));
#use fzero to get Berry-Esseen tolerance
tol1 = min(NCheb_inv,NBE_inv);
#take the min of Chebyshev and Berry Esseen tolerance
return(tol1)
}

meanMC_g_param = function(Yrand,abstol,reltol,nSig,n1,alpha,fudge,tbudget,nbudget) {
  out_param.Yrand = Yrand
  out_param.abstol = abstol
  out_param.reltol = reltol
  out_param.nSig = nSig
  out_param.n1 = n1
  out_param.alpha = alpha
  out_param.fudge = fudge
  out_param.tbudget = tbudget
  out_param.nbudget = nbudget
  out_param = c(out_param.Yrand,out_param.abstol,out_param.reltol,out_param.nSig,out_param.n1,out_param.alpha,out_param.fudge,out_param.tbudget,out_param.nbudget)
  return(out_param)
}

meanMC_g_err = function(out_param)
  # Handles errors in meanMC_g and meanMC_g_param to give an exit with
  #  information.
  #            out_param.exit = 0   success
  #                             1   too many samples required
  
  if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
case 1 % not enough samples to estimate the mean.
nexceed = out_param.n(out_param.tau);
warning('MATLAB:meanMC_g:maxreached',...
        [' In order to achieve the guaranteed accuracy, at step '...
         int2str(out_param.tau) ', tried to evaluate at ' int2str(nexceed) ...
         ' samples, which is more than the remaining '...
         int2str(out_param.nremain) ...
         ' samples. We will use all the samples left to estimate the mean without guarantee.']);
return
  }
}
