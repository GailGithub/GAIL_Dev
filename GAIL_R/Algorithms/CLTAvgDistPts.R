# Average Distance Between Points Using CLT Confidence Intervals

# A simple, concrete example  
# Recall again the problem of determining the average distance between two
# points in a unit square.  The true answer is a 4-dimensional integral
# 
# \[ \mu = \int_{[0,1]^2 \times [0,1]^2} \sqrt{(x_1-y_1)^2 +
# (x_2-y_2)^2} \, {\rm d}x_1 \, {\rm d}x_2 \, {\rm d}y_1 \, {\rm d}y_2\]
#                                             
# We talked about how to compute the answer using Monte Carlo methods.
# First we define a mean distance function as follows

distfun = function(n){sqrt(rowSums((matrix(runif(n*2),n,2) - matrix(runif(n*2),n,2))^2))}

# 
# To approximate the answer, \(\mu\), we may use the sample mean with |n =
# 1e6| points.  This is the Monte Carlo method
tstart = proc.time()
meandist = mean(distfun(1e6))
proc.time() - tstart

#
# The problem is that we do not yet know how accurate our answer is.

# Central Limit Theorem (CLT) confidence intervals
# The probability distribution of the sample mean with \(n\) samples of IID
# \(Y_i\) approaches the Gaussian (or normal) distribution with mean $E(Y)$
# and variance $\mbox{var}(Y)/n$.  This idea is behind the algorithm
# |meanMC_CLT|

# Need to store meanMC_CLT
source("meanMC_CLT.R")
# 
# We may use this algorithm to find a fixed-width confidence interval
tstart1 = proc.time()
meandist = meanMC_CLT(distfun,0.02)
proc.time() - tstart1

#
# Note that the time required by |meanMC_CLT| is less than that required by
# taking the mean of 1e6 samples.  To find out how many samples
# |meanMC_CLT| actually used, we may try

meandist = meanMC_CLT(distfun,0.02)

#
# The value of |output.ntot| is the total number of samples used.

# Sometimes the CLT confidence interval does not work
# Unfortunately, |meanMC_CLT| does not have solid theoretical support.
# E.g. consider the example

a=1e2; # a parameter
Y = function(n){rnorm(n) + a*(a*runif(n)<1)} #a mixture distribution with mean 1
for (i in 1:4){
muhat = meanMC_CLT(Y,0.01)
print(muhat)} 
#try out multiple times with tolerance 0.01
#
# Note that the answers \(\pm\) 0.01 do not overlap.  Thus, they must be
# wrong. The problem is that the number of samples used to estimate the
# variance (default = 100) is too small for this distribution with a large
# kurtosis. If we increase the number of samples used to estimate the
# variance to 10000, then the answers are correct.

for (i in 1:4) {
muhat = meanMC_CLT(Y,0.01,0.01,10000) #try out multiple times with tolerance 0.01
print(muhat)
}

#
# The true mean in this case is 1.
#
# This is great, but how does one know that 10000 samples for estimating
# the variance is enough, but 100 is too few.  We would like a theorem.  If
# we want to have algorithms with solid theoretical justification, then we
# should look at the *Guaranteed Automatic Integration Library (GAIL)*.

# Authors: Anthony Karahalios and Luana Terra