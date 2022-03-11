a = 1/sqrt(2);  
d = 3;
abstol = 0.005;
reltol = 0;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting Inf
yinv = @(t)(erfcinv(replaceZeros(abs(t))));
f1 = @(t,d) cos(sqrt(normsqd(yinv(t))))*(sqrt(pi))^d;
fKeister = @(x) f1(x,d); 
inputArgs = {'absTol',abstol,'relTol',reltol};
hyperbox = [zeros(1,d); ones(1,d)];
[u1,~] = cubMC_g(fKeister, hyperbox, inputArgs{:});
[u2,~] = cubSobol_g(fKeister, hyperbox, inputArgs{:});
[u3,~] = cubLattice_g(fKeister, hyperbox, inputArgs{:});
[~,u4] = cubBayesNet_g(fKeister, d, inputArgs{:});
[~,u5] = cubBayesLattice_g(fKeister, d, inputArgs{:});