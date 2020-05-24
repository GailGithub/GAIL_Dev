%DT_CUBBAYESLATTICE_G  lengthy doctests for cubBayesLattice_g
%
%
%   Example 1: Multivariate normal probability
%   Steepest gradient descent used with analytic kernel gradient in
%   parameter search
% 
%   >> dim=2; absTol=1e-3; relTol=1e-2; fName = 'MVN';
%   >> C = [4 1 1; 0 1 0.5; 0 0 0.25]; MVNParams.Cov = C'*C; MVNParams.C = C;
%   >> MVNParams.a = [-6 -2 -2]; MVNParams.b = [5 2 1]; MVNParams.mu = 0;
%   >> MVNParams.CovProp.C = chol(MVNParams.Cov)';
%   >> muBest = 0.676337324357787;
%   >> integrand =@(t) GenzFunc(t,MVNParams);
%   >> inputArgs={'absTol',absTol,'relTol',relTol};
%   >> inputArgs=[inputArgs {'order',1,'ptransform','C1sin','arbMean',true}];
%   >> inputArgs=[inputArgs {'useGradient',true}];
%   >> [obj,muhat]=cubBayesLattice_g(integrand,dim, inputArgs{:});
%   >> check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))
%   check = 1
% 
% 
%   Example 2: Keister function
%   Kernel order r chosen automatically
% 
%   >> dim=2; absTol=1e-3; relTol=1e-2;
%   >> normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
%   >> replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
%   >> yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
%   >> ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
%   >> fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
%   >> inputArgs ={'absTol',absTol, 'relTol',relTol};
%   >> inputArgs =[inputArgs {'order',0, 'ptransform','C1','arbMean',true}];
%   >> obj=cubBayesLattice_g(fKeister,dim,inputArgs{:});
%   >> [muhat,outParams] = compInteg(obj);
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%   >> check = double(outParams.optParams.r > 0)
%   check = 1
% 
% 
%   Example 3: Another example using dimension specific shape parameter
%
%   >> const = [1E-4 1 1E4];
%   >> fun = @(x)sum(bsxfun(@times, const, sin(2*pi*x.^2)), 2);
%   >> dim=3; absTol=1e-3; relTol=1e-2;
%   >> exactInteg = fresnels(2)*sum(const)/2;
%   >> inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
%   >> inputArgs = [{'absTol',absTol,'oneTheta',false,'useGradient',false} inputArgs];
%   >> obj=cubBayesLattice_g(fun, dim, inputArgs{:});
%   >> [muhat,outParams]=compInteg(obj);
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%   >> etaDim = size(outParams.optParams.aMLEAll, 2)
%   etaDim = 3
%
%
%   Example 4: Legacy API with only two arguments: integrand and number of dimensions
%   Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over
%   the interval [0,1] with parameters: order=1, abstol=0.001, relTol=0.01
% 
%   >> fun = @(x) exp(sum(cos(2*pi*x), 2));
%   >> dim=2; 
%   >> exactInteg = besseli(0,1)^dim;
%   >> [obj,muhat]=cubBayesLattice_g(fun,dim);
%   >> check = double(abs(exactInteg-muhat) < 0.01)
%   check = 1
%
%
%   Example 5: Keister function
%   Same shape parameter for all dimensions
% 
%   >> dim=3; absTol=1e-3; relTol=1e-2;
%   >> normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
%   >> replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
%   >> yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
%   >> ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
%   >> fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
%   >> inputArgs ={'absTol',absTol, 'relTol',relTol};
%   >> inputArgs =[inputArgs {'order',2, 'ptransform','C1','arbMean',true}];
%   >> obj=cubBayesLattice_g(fKeister,dim,inputArgs{:});
%   >> [muhat,outParams]=compInteg(obj);
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%   >> etaDim = size(outParams.optParams.aMLEAll, 2)
%   etaDim = 1
%
%
%   Example 6:
%   Estimate the integral with integrand Genz "Oscillatory" in
%   genz_test_fun, index 1.
%   >> in_param.dim = 3;alpha = ones(1,in_param.dim); beta = 1./ (1:in_param.dim);r=1;
%   >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 1;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.06***
%   >> [obj,Q] = cubBayesLattice_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-2);
%   >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-2,1,f_true,'max'))
%   check = 1
%  
%  
%   Example 7:
%   Estimate the integral with integrand Genz "Product Peak" in
%   genz_test_fun, index 2.
%   >> index = 2;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.66***
%   >> [~,Q] = cubBayesLattice_g(f,in_param.dim, 'absTol',1e-3, 'relTol',1e-2);
%   >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-3,1,f_true,'max'))
%   check = 1
%  
%  
%   Example 8:
%   Estimate the integral with integrand Genz "Product Peak" in
%   genz_test_fun, index 3.
%   >> index = 3;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.04***
%   >> [~,Q] = cubBayesLattice_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-2);
%   >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-2,1,f_true,'max'))
%   check = 1
%  
%
%   Example 9:
%   Estimate the integral with integrand Genz "Gaussian" in
%   genz_test_fun, index 4.
%   >> index = 4;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.62***
%   >> [~,Q] = cubBayesLattice_g(f,in_param.dim, 'absTol',1e-3, 'relTol',1e-3);
%   >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-3,1,f_true,'max'))
%   check = 1
%  
%  
%   Example 10:
%   Estimate the integral with integrand Genz "Continuous" in
%   genz_test_fun, index 5.
%   >> index = 5;
%   >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
%   >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
%   f_true = 0.38***
%   >> [~,Q] = cubBayesLattice_g(f, in_param.dim, 'absTol',1e-3, 'relTol',1e-3);
%   >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-3,1,f_true,'max'))
%   check = 1
%    
% 
%   Example 11: 
%   Demonstrates passing the optional params as a struct
% 
%   Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over the
%   interval [0,1]^2 with parameters: order=2, ptransform=C1sin, abstol=0.001,
%   relTol=0.01
% 
%   >> fun = @(x) exp(sum(cos(2*pi*x), 2));
%   >> dim=2; absTol=1e-2; relTol=1e-2;
%   >> exactInteg = besseli(0,1)^dim;
%   >> inParams.absTol = absTol;
%   >> inParams.relTol = relTol;
%   >> inParams.oneTheta = false;
%   >> inParams.ptransform = 'C1sin';
%   >> inParams.order = 2;
%   >> obj=cubBayesLattice_g(fun, dim, inParams);
%   >> [muhat,outParams]=compInteg(obj);
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%   >> etaDim = size(outParams.optParams.aMLEAll, 2)
%   etaDim = 2
%
%
%   Example 12: 
%   Demonstrates 4 arguments passing style
% 
%   Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over the
%   interval [0,1]^2 with tolerance: abstol=0.001, relTol=0.01
% 
%   >> fun = @(x) exp(sum(cos(2*pi*x), 2));
%   >> dim=2; absTol=1e-2; relTol=1e-2;
%   >> exactInteg = besseli(0,1)^dim;
%   >> obj=cubBayesLattice_g(fun, dim, absTol, relTol);
%   >> [muhat,outParams]=compInteg(obj);
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%   >> etaDim = size(outParams.optParams.aMLEAll, 2)
%   etaDim = 1
 