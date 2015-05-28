stdnormcdf = function (z) {
  # this function is to define cumulative distribution function (CDF) of the
  # standard normal distribution.
  p = 0.5*(2*pnorm(-z/sqrt(2)*sqrt(2),lower=FALSE))
  # Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
  # to produce accurate near-zero results for large negative x.
  return(p)
}
