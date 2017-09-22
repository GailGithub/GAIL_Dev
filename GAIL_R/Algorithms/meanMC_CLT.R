#MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
#
#   tmu = MEANMC_CLT(Yrand,abstol,alpha,nSig,fudge) estimates the mean, mu, of a random variable Y to
#   within a specified error tolerance, i.e., | mu - tmu | <= abstol with
#   probability at least 1-alpha, where abstol is the absolute error
#   tolerance.  The default values are abstol=0.01 and alpha=0.01. Input
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
#     function(n) {runif(n)}.
#
#     abstol --- the absolute error tolerance, which should be
#     positive, default value is 0.01.
#
#     alpha --- the uncertainty, which should be a small positive
#     percentage. The default value is 0.01.
#
#     nSig --- the number of samples used to compute the sample variance
#
#     fudge --- the standard deviation invlation factor
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
#Default Values
#fudge = 1.2; variance inflation factor
#nSig = 1e2; number of samples to estimate variance  100
#alpha = 0.01; uncertainty  
#abstol = 0.01; absolute error tolerance  
#Yrand = @(n) rand(n,1); random number generator 
#
# Authors: Anthony Karahalios, Luana Terra, Marcela Ribeiro, Ramon Oliveira, Joao Mateus Cunha
# Youhan Lu, Batsuren Davaademberel

meanMC_CLT = function(Yrand = function(n) {runif(n)},abstol = 0.01,alpha = 0.01,nSig = 1e2,fudge = 1.2) {
  
check=meanMC_g_param(Yrand,abstol,alpha,nSig,fudge) 
  
nMax=1e8; #maximum number of samples allowed.
Yrand=check$Yrand;
alpha = check$alpha; #save the input parameters to a structure
fudge = check$fudge;
nSig = check$nSig;
abstol = check$abstol;
start <- Sys.time (); #start the clock 
Yval = Yrand(nSig);# get samples to estimate variance 
out_param.var = var(Yval); #calculate the sample variance--stage 1
sig0 = sqrt(out_param.var); #standard deviation
sig0up = fudge*sig0; #upper bound on the standard deviation
hmu0 = mean(Yval); # FOI ADICIONADO
nmu = max(1,ceiling((-qnorm(alpha/2)*sig0up/max(abstol,alpha*abs(hmu0)))^2)); # Changed
#number of samples needed for mean
stopifnot(nmu<nMax) 
#don't exceed sample budget
tmu = mean(Yrand(nmu)); #estimated mean
out_param.ntot = nSig + nmu; #total samples required
out_param.time = Sys.time () - start;#proc.time() - tstart; #elapsed time
#out_param.time = unname(out_param.time)
output = c("tmu" = tmu,"out_param.ntot" = out_param.ntot,"out_param.var" = out_param.var,"out_param.time" = out_param.time)
return(output)
}

meanMC_g_param = function(Yrand,abstol,alpha,nSig,fudge) {
  
  
  if(!is.function(Yrand)) {message("Yrand must be a function - Now meanMC_CLT is using default Yrand=function(n) {runif(n)}")
                           Yrand=function(n) {runif(n)}
  }
  if(abstol<=0) {message("Absolute error tolerance should be greater than 0 - Now meanMC_CLT is using default absolute error tolerance = 0.01.")
                 abstol = 0.01
  }
  if(alpha<=0 | alpha>=1) {message("The uncertainty should be between 0 and 1 - Now meanMC_CLT is using default value = 0.01.")
                           alpha = 0.01
  }
  if(nSig%%1!=0 | nSig<30) {message("The number nSig should be a positive integer at least 30 - Now meanMC_CLT is using default value = 1e2.")
                            nSig = 1e2
  }
  if(fudge <= 1) {message("'The fudge factor should be larger than 1 - Now meanMC_CLT is using default value = 1.2.")
                  fudge = 1.2
  }
  out_param=c("Yrand"=Yrand,"abstol"=abstol,"nSig"=nSig,"alpha"=alpha,"fudge"=fudge)
  return(out_param)                      
  
}
