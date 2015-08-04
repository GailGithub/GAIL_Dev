tolfun = function(abstol,reltol,theta,mu,toltype) {
  # TOLFUN generalized error tolerance function.
  #
  # Input Parameters:
  # abstol --- absolute error tolertance
  # reltol --- relative error tolerance
  # theta --- parameter in 'theta' case
  # mu --- true mean
  # toltype --- different options of tolerance function
  
  if (toltype == "comb") {    # the linear combination of two tolerances
    tol  = theta*abstol+(1-theta)*reltol*abs(mu)
  }
  else if (toltype == "max") {    # the max case
    tol  = max(abstol,reltol*abs(mu))
  }
  #theta=0---relative error tolarance
  #theta=1---absolute error tolerance
  return(tol)
}