test1.in_param.abstol = 2e-4;
test1.in_param.reltol = 0;
test1.Yrand = @(n) exp(rand(n, 1));


k1 = meanMC(test1);

tic,genMu(k1),toc


test2.in_param.abstol = 2e-4;
test2.in_param.reltol = 0;
test2.YXrand = @expr;
test2.muX = 0.5;
test2.method = 'cv';

k2 = meanMC(test2);

tic,genMu(k2),toc