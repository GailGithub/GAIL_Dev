test1.in_param.abstol = 2e-4;
test1.in_param.reltol = 0;
test1.Yrand = @(n) exp(rand(n, 1));


k1 = meanMC(test1);

tic,genMu(k1),toc


test2.in_param.abstol = 2e-4;
test2.in_param.reltol = 0;
test2.cv_param.YXrand = @expr;
test2.cv_param.muX = 0.5;
test2.cv_nb = 1e3;
test2.method = {'cv'};

k2 = meanMC(test2);

[a,b] = genMu(k2);


distfun1 = @distfun;
distfun2 = @(n) sqrt(sum((rand(n,2)  - rand(n,2)).^2,2));

test3.in_param.abstol = 2e-4;
test3.in_param.reltol = 0;
test3.cv_param.YXrand = distfun1;
test3.cv_param.muX = repmat(0.5,1,4);
test3.Yrand = distfun2;
test3.method = {'cv','plain'};
k3 = meanMC(test3);

tic,genMu(k3),toc
