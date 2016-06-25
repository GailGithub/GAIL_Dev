obj1.in_param.abstol = [2e-4 1];
obj1.in_param.reltol = 0;
obj1.method = {'plain'}; 
obj1.Yrand = {@(n) sin(rand(n, 1)), @(n) exp(rand(n, 1))};
test1 = meanMC(obj1);
genMu(test1)