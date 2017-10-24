#MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
#
#   sol = MEANMC_CLT(Y,exact,absTol,relTol,alpha,nSig,inflate) estimates the
#   mean, mu, of a random variable to within a specified error tolerance,
#   i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
#   1-alpha, where abstol is the absolute error tolerance.  The default
#   values are abstol=1e-2 and alpha=1#. Input Y is a function handle that
#   accepts a positive integer input n and returns an n x 1 vector of IID
#   instances of the random variable.
#
#   This is a heuristic algorithm based on a Central Limit Theorem
#   approximation.
#
#
#   Input Arguments
#
#     Y --- the function or structure for generating n IID instances of a
#     random variable Y whose mean we want to estimate. Y is often defined
#     as a function of some random variable X with a simple distribution.
#     The input of Yrand should be the number of random variables n, the
#     output of Yrand should be n function values. For example, if Y = X.^2
#     where X is a standard uniform random variable, then one may define
#     Yrand = @(n) rand(n,1).^2.
#     
#
#     absTol --- the absolute error tolerance, which should be
#     non-negative --- default = 1e-2
#
#     relTol --- the relative error tolerance, which should be
#     non-negative and no greater than 1 --- default = 0
#
#     alpha --- the uncertainty, which should be a small positive
#     percentage --- default = 1#
  #
#     nSig --- the number of samples used to compute the sample variance
#     --- default = 1000
#
#     inflate --- the standard deviation inflation factor --- default = 1.2
#
#   Output Arguments
#
#     Y --- the random generator
#
#     absTol --- the absolute error tolerance
#
#     relTol --- the relative error tolerance
#
#     alpha --- the uncertainty
#
#     mu --- the estimated mean of Y.
#
#     stddev --- sample standard deviation of the random variable
#
#     nSample --- total sample used.
#
#     time --- the time elapsed in seconds.
#
#     errBd --- the error bound.
#
# Youhan Lu, Batsuren Davaademberel
#This is the example 1
#Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2 with absolute tolerance 1e-3 and relative tolerance 0:
library(pracma)
# YXn<-(function(n) (rand(n,1)^2))
# nY=1
# true_value=list()
# exact<-1/3
# s<-list(YXn,c(nY),c(true_value))
# names(s)<-c('Y','nY','trueMuCV')

#This is the example 2
#Estimate the integral f(x)=exp(-x^2) in the interval [0,1] using x as a control variate and relative error 1e-3:
# f<-(function(x)(list(exp(-x^2),x)))
# YXn<-(function(n) (f(rand(n,1))))
# nY<-1
# true_value<-list(1/2)
# exact<-erf(1)*sqrt(pi)/2
# s<-list(YXn,c(nY),c(true_value))
# names(s)<-c('Y','nY','trueMuCV')
# meanMC_CLT(s,exact)


#This is the example 3
#Estimate the Keister's integration in dimension 1 with a=1, 1/sqrt(2) and using cos(x) as a control variate:
normsqd = (function(x) (rowSums(x*x)));
f<-(function(normt,a,d) ((2*pi*a^2)^(d/2)) * cos(a*sqrt(normt))* exp((1/2-a^2)*normt))
f1<-(function(x,a,d) f(normsqd(x),a,d))
f2<-(function(x) (list((f1(x,1,1)),(f1(x,1/sqrt(2),1)),cos(x))))
YXn<-(function(n) (f2(randn(n,1))))
true_value<-list(1/sqrt(exp(1)))
nY<-2
exact<-1.380388447043143
s<-list(YXn,c(nY),c(true_value))
names(s)<-c('Y','nY','trueMuCV')
meanMC_CLT(s,exact)


# This is the example 4
# f<-(function(x) (list((x[,1]^3*x[,2]^3*x[,3]^3),(x[,1]*x[,2]*x[,3]))))
# YXn<-(function(n) (f(rand(n,3))))
# true_value<-list(1/8)
# nY<-1
# exact<-1/64
# s<-list(YXn,c(nY),c(true_value))
# names(s)<-c('Y','nY','trueMuCV')
# meanMC_CLT(s,exact)


# This is the example 5
# f<-(function(x) (list((x[,1]^3*x[,2]^3*x[,3]^3),(x[,1]^2*x[,2]^2*x[,3]^2-1/27+1/64),(x[,1]*x[,2]*x[,3]),(x[,1]+x[,2]+x[,3]))))
# YXn<-(function(n) (f(rand(n,3))))
# true_value<-list(1/8,1.5)
# nY<-2
# s<-list(YXn,c(nY),c(true_value))
# names(s)<-c('Y','nY','trueMuCV')
# meanMC_CLT(s,exact)
# This is the main function 
meanMC_CLT = function(s,exact,abstol = 1e-3,relTol=0.0, alpha = 1e-2,nSig = 1e2,inflate = 1.2) 
{
  
  nMax=16777216
  pmt<-proc.time() # start the clock, pmt is the current time
  Yrand<-s$Y # the random number generator
  q<-s$nY    # the number of target random variable
  true_value<-s$trueMuCV
  p<-length(true_value) # the number of control variates
  print(c(p,q))
  val<-Yrand(nSig) # get samples to estimate variance
if (p==0 && q==1){
    YY<-val
} else {
    # if there is control variate, construct a new random variable that has the same expected value and smaller variance
    meanVal<-t(matrix(colMeans(data.frame(val)))) # the mean of each column 
    val_matrix<-do.call(cbind,val) # convert val into matrix thus we can perform matrix transformation 
    A<-sweep(val_matrix, 2, meanVal, '-') # performing Value minus mean row by row
    svdResult<-svd(rbind(A,t(c(ones(1,q),zeros(1,p)))),LINPACK = FALSE) #use SVD to solve a constrained least square problem
    Sdiag<-diag(svdResult$d) # the vector of the single values
    U2<-t(svdResult$u[dim(svdResult$u)[1],]) # the last row of U 
    beta<-svdResult$v%*%inv(Sdiag)%*%t(U2)%*%(inv(U2%*%t(U2))) # get the coefficient for control variates
    YY<-data.matrix(data.frame(data.frame(val)[1:q],data.frame(A)[(q+1):dim(A)[2]]))%*% beta #get samples of the new random variable 
}
  stddev<-std(YY) # standard deviation of the new samples
  sig0up<-inflate*stddev #upper bound on the std
  hmu0<-mean(YY) # mean of the samples
  nmu=max(1,ceil((sqrt(2)*erfcinv(2*alpha/2)/max(abstol*abs(hmu0)))^2))
  #number of samples needed for the error tolerance
  if (nmu>nMax) {# exceed sample budget
    warning(capture.output(cat('The algorithm wants to use nmu=',nmu,', which is too big. Using ',nMax,' instead.')))
    nmu=nMax # revise nmu
  }
  print(p)
  YY=Yrand(nmu) #get samples for computing the mean
  if (p>0 ||q>1)
  {
    if (length(true_value)!=0){
      print('ture_value exists')
      YY_matrix<-do.call(cbind,YY) # convert YY into matrix thus we can perform matrix transformation 
      YY<-YY_matrix
      true_value_matrix<-do.call(cbind,true_value) # matrix transformation from list
      if ((q+1)==dim(YY_matrix)[2]){
        YY[,(q+1)]<-sweep(matrix(YY_matrix[,(q+1)],ncol=1),2,true_value_matrix,'-')
      }else{
      YY[,((q+1):dim(YY_matrix)[2])]<-sweep(matrix(YY_matrix[,((q+1):dim(YY_matrix)[2])],ncol=dim(YY_matrix)[2]), 2, true_value_matrix, '-') # performing Value minus mean row by row
      }
    }
  YY<-YY%*%beta #incorporate the control variates and multiple Y's
  }
  sol<-colMeans(YY) #estimated mean
  nSample<-nSig+nmu #total samples required
  time<-proc.time()-pmt #total time for calculating MC
  errBd<-(sqrt(2)*erfcinv(2*alpha/2))*sig0up/sqrt(nmu)
  print(hmu0)
  check<-(abs(exact-sol)<max(0,1e-3*abs(exact)))
  print(check)
  return (sol)
}
sol<-meanMC_CLT(s,exact)
print(exact)
print(sol)
