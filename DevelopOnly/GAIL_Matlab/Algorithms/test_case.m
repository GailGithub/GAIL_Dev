tic, [tmu,out_param] = meanMCCV_g(@distfun, [0.5 0.5 0.5 0.5], 0.0002, 0), toc

distfun2 = @(n) sqrt(sum((rand(n,2)  - rand(n,2)).^2,2));
tic, [tmu,out_param] = meanMC_g(distfun2,0.0002,0), toc
