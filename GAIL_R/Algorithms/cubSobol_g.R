####sobstr=sobolset(out_param.d); %generate a Sobol' sequence
####sobstr=scramble(sobstr,'MatousekAffineOwen'); %scramble it
#sobstr=sobol(n,out_param.d,scrambling=1)
Stilde=rep(0,out_param.mmax-out_param.mmin+1); #initialize sum of DFWT terms
#CStilde_low = -inf(1,out_param.mmax-l_star+1); #initialize #various sums of DFWT terms for necessary conditions
Cstilde_low=matrix(-Inf,1,out_param.mmax-l_star+1);
#CStilde_up = inf(1,out_param.mmax-l_star+1); #initialize various #sums of DFWT terms for necessary conditions
CStilde_up = matrix(Inf,1,out_param.mmax-l_star+1);
errest=rep(0,out_param.mmax-out_param.mmin+1); #initialize error estimates
appxinteg=rep(0,out_param.mmax-out_param.mmin+1); #initialize approximations to integral
exit_len = 2;
#out_param.exit=false(1,exit_len); #we start the algorithm with #all warning flags down