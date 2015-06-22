library(pracma)
#MEANMCabs_g_abbr Monte Carlo method to estimate the mean of a random variable
#
#   tmu = MEANMCabs_g_abbr(Yrand,abstol,alpha,nSig,fudge) estimates the mean, mu, of a random variable Y to
#   within a specified error tolerance, i.e., | mu - tmu | <= abstol with
#   probability at least 1-alpha, where abstol is the absolute error
#   tolerance.  The default values are abstol=1e-2 and alpha=1%. Input
#   Yrand is a function handle that accepts a positive integer input n and
#   returns an n x 1 vector of IID instances of the random variable Y.
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
#     abstol --- the absolute error tolerance, which should be
#     positive, default value is 1e-2.
#
#     alpha --- the uncertainty, which should be a small positive
#     percentage. The default value is 1%.
#
#     nSig --- the number of samples used to compute the sample variance. The
#     default value is 1e2
#
#     fudge --- the standard deviation invlation factor.  The default value
#     is 1e2
#
#   Output Arguments
#
#     tmu --- the estimated mean of Y.
#
#     out_param.ntot --- total sample used.
#
#     out_param.var --- the sample variance.
#
#     out_param.time --- the time elapsed in seconds.
#

#This is a heuristic algorithm based on a Central Limit Theorem
#approximation

meanMCabs_g_abbr = function (Yrand = function(n) {runif(n)^2},abstol=0.01,alpha=0.01,nSig=100,fudge=1.2){
  
  check=meanMC_g_param(Yrand,abstol,alpha,nSig,fudge)  
  
  nMax=1e8; #maximum number of samples allowed.
  Yrand=check$Yrand;
  out_param.alpha = check$alpha; #save the input parameters to a structure
  out_param.fudge = check$fudge;
  out_param.nSig = check$nSig;
  out_param.abstol = check$abstol;
  tstart = proc.time(); #start the clock
  Yval = Yrand(nSig); # get samples to estimate variance
  out_param.var = var(Yval); #calculate the sample variance--stage 1
  sig0 = sqrt(out_param.var); #standard deviation
  sig0up = out_param.fudge*sig0; #upper bound on the standard deviation
  alpha1 = 1-sqrt(1-out_param.alpha); # one side of the uncertainty
  out_param.kurtmax = (out_param.nSig-3)/(out_param.nSig-1) 
  + ((alpha1*out_param.nSig)/(1-alpha1))*(1-1/out_param.fudge^2)^2;
  toloversig = out_param.abstol/sig0up;
  # absolute error tolerance over sigma
  ncheb = ceiling(1/(alpha1*toloversig^2));
  # use Chebyshev inequality to estimate n
  A=18.1139;
  A1=0.3328;
  A2=0.429; # three constants in Berry-Esseen inequality
  M3upper=out_param.kurtmax^(3/4); #using Jensen inequality to 
  #                                       bound the third moment
  
  BEfun= function(logsqrtn) {pnorm(-exp(logsqrtn)*toloversig)+
                               exp(-logsqrtn)*min(A1*(M3upper+A2),
                                                  A*M3upper/(1+(exp(logsqrtn)*toloversig)^3))- alpha1/2;
  }
  # Berry-Esseen Inequality
  logsqrtnCLT=log(qnorm(1-alpha1/2)/toloversig);
  nmu=min(ncheb,ceiling(exp(2*fzero(BEfun,logsqrtnCLT)$x)));
  #get the min n (used to estimate mu) by using cheb and BEfun
  if(nmu>=nMax) {message("Program stopped due to the guarantee needing over 1e8 samples")}
  stopifnot(nmu<nMax) 
  #don't exceed sample budget
  tmu=mean(Yrand(nmu)); #estimated mean
  out_param.ntot=nSig+nmu; #total samples required
  out_param.time=proc.time()-tstart; #elapsed time
  out_param.time=unname(out_param.time)
  return(c("tmu"= tmu,"out_param.ntot" = out_param.ntot,"out_param.var" = out_param.var,"out_param.time" = out_param.time[3]))
}


meanMC_g_param = function(Yrand,abstol,alpha,nSig,fudge) {
  
  
  if(!is.function(Yrand)) {message("Yrand must be a function - Now R is using default Yrand=function(n) {runif(n)^2}")
                           Yrand=function(n) {runif(n)^2}
  }
  if(abstol<=0) {message("Absolute error tolerance should be greater than 0 - Now R is using default absolute error tolerance = 1e-2.")
                 abstol = 1e-2
  }
  if(nSig%%1!=0 | nSig<30) {message("The number nSig should be a positive integer at least 30 - We will use the default value 1e2.")
                                            nSig = 1e2
  }
  if(alpha<=0 | alpha>=1) {message("The uncertainty should be between 0 and 1 - We will use the default value 0.01.")
                           alpha = 0.01
  }
  if(fudge <= 1) {message("'The fudge factor should be larger than 1 - We will use the default value 1.2.")
                  fudge = 1.2
  }
  out_param=c("Yrand"=Yrand,"abstol"=abstol,"nSig"=nSig,"alpha"=alpha,"fudge"=fudge)
  return(out_param)                      
  
}
