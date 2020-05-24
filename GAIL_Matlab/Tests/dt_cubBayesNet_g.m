%DT_CUBBAYESNET_G  lengthy doctests for cubBayesNet_g
%
%
%   Example 1: 
%   Estimate the integral with integrand Genz "Oscillatory" in
%   genz_test_fun, index 1.
%   >>  in_param.dim = 3;alpha = ones(1,in_param.dim); beta = 1./ (1:in_param.dim);r=1;
%   >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 1;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
%   f_true = 0.06***
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-2); Q = Q + 0.001
%   Q = 0.06***
% 
%
%   Example 2: 
%   Estimate the integral with integrand Genz "Product Peak" in
%   genz_test_fun, index 2.
%   >> index = 2;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
%   f_true = 0.66***
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-3); Q = Q + 0.001
%   Q = 0.66***
%   
%
%   Example 3: 
%   Estimate the integral with integrand Genz "Product Peak" in
%   genz_test_fun, index 3.
%   >> index = 3;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
%   f_true = 0.04***
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 1e-3, 1e-2); Q = Q + 0.001
%   Q = 0.04***
%   
%
%   Example 4: 
%   Estimate the integral with integrand Genz "Gaussian" in
%   genz_test_fun, index 4.
%   >> index = 4;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
%   f_true = 0.62***
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 1e-3, 1e-3); Q = Q + 0.001
%   Q = 0.62***
%  
%
%   Example 5: 
%   Estimate the integral with integrand Genz "Continuous" in
%   genz_test_fun, index 5.
%   >> index = 5;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.38***
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-3)
%   Q = 0.38***
%  
%
%   Example 6: 
%   Estimate the integral with integrand Genz "Discontinuous" in
%   genz_test_fun, index 6.
%   >> index = 6;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.44***
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-3)
%   Q = 0.44***
%  
%
%   Example 7: 
%   Estimate the integral with integrand Keister test function in
%   genz_test_fun, index 7.
%   >> index = 7; in_param.dim = 3;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> [~,Q] = cubBayesNet_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-4)
%   Q = 2.16***
%  
%  
%   Example 8: Multivariate normal probability
%   >> dim=2; absTol=1e-3; relTol=1e-2; fName = 'MVN';
%   >> C = [4 1 1; 0 1 0.5; 0 0 0.25]; MVNParams.Cov = C'*C; MVNParams.C = C;
%   >> MVNParams.a = [-6 -2 -2]; MVNParams.b = [5 2 1]; MVNParams.mu = 0;
%   >> MVNParams.CovProp.C = chol(MVNParams.Cov)';
%   >> muBest = 0.676337324357787;
%   >> integrand =@(t) GenzFunc(t,MVNParams);
%   >> inputArgs={'absTol',absTol,'relTol',relTol, 'arbMean',true};
%   >> obj=cubBayesNet_g(integrand,dim, inputArgs{:});
%   >> [muhat,outParams] = compInteg(obj);
%   >> check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))
%   check = 1
%   >> etaDim = size(outParams.optParams.aMLEAll, 2)
%   etaDim = 1
%
%
%   Example 9: Legacy API with only two arguments: integrand and number of dimensions
%   Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over
%   the interval [0,1] with parameters: order=1, abstol=0.001, relTol=0.01
% 
%   >> fun = @(x) exp(sum(cos(2*pi*x), 2));
%   >> dim=2; 
%   >> exactInteg = besseli(0,1)^dim;
%   >> [~,muhat]=cubBayesNet_g(fun,dim);
%   >> check = double(abs(exactInteg-muhat) < 0.01)
%   check = 1
% 
% 
%   Example 10: Keister function
% 
%   >> dim=2; absTol=1e-3; relTol=1e-2;
%   >> normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
%   >> replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
%   >> yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
%   >> ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
%   >> fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
%   >> inputArgs ={'absTol',absTol, 'relTol',relTol};
%   >> inputArgs =[inputArgs {'arbMean',true}];
%   >> [obj,muhat]=cubBayesNet_g(fKeister,dim,inputArgs{:});
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
% 
% 
