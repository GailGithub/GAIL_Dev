% % Example 1:
% % Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:
% >> w.f = @(x) prod(x,2);
% >> w.absTol=1e-3;
% >> w.relTol=0;
% >> w.relTol=0; 
% >> w.domain = [zeros(1,2);ones(1,2)];
% >> [q, out] = cubMC_CLT(w)
% exactsol = 1/4;
% 
% % Example 2:
% >> w.f= @(x) exp(-x(:,1).^2-x(:,2).^2); w.absTol=1e-3;
% >> w.relTol=1e-2; w.domain = [-ones(1,2);2*ones(1,2)];
% >> [q, out_param] = cubMC_CLT(w);
% exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
% 
% % Example 3: 
% >> w.f= @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2]; w.domain = [zeros(1,3);2*ones(1,3)];
% >> w.trueMuCV=[8,32/3]; w.absTol=1e-3;
% >> w.relTol=0; w.measure = 'uniform';
% >> [q, out_param] = cubMC_CLT(w);
% exactsol = 128/3;
% 
% % Example 4:
% >> w.f= @(x) 3./(5-4*(cos(2*pi*x)));
% >> w.absTol=1e-3; w.relTol=0;
% >> w.domain = [0;1]; 
% >> [q, out_param] = cubMC_CLT(w);
% exactsol = 1;
% 
% % Example 5: 
% >> w.f= @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2]; w.domain = [zeros(1,3);2*ones(1,3)];
% >> w.trueMuCV=[8,32/3]; w.absTol=1e-3;
% >> w.relTol=0; w.measure = 'uniform';
% >> [q, out_param] = cubMC_CLT(w);
% exactsol = 128/3;
