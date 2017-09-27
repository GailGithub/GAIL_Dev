function [integ, out] = cubMC_CLT(varargin)
% CUBMC_CLT approximates integrals via Monte Carlo simulation with Central
% Limit Theorem (CLT) confidence intervals constructed to fit the error
% tolerance
%
%   [integ, out] = CUBMC_CLT(f, domain, domainType, measure)
% 
% The main work is done by MEANMC_CLT.  Input can also be done via
% structures or classes.
%
% Input Arguments
%   f --- the integrand, in general a function of d variables
%
%   domain --- the integration domain represented as a 2 x d array
%              if a box, the top row is the lower left corner and the
%              bottom row is the top right corner; if a ball or sphere, the
%              top row is the center and the first element of the bottom
%              row is the radius
%               
%   domainType --- 'box', 'ball', or 'sphere'
%   
%   measure --- integration measure: 'uniform', 'Gaussian', or 'Lebesgue'
%
% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2
% >> integ = cubMC_CLT(@(x) prod(x,2), [0 0; 1 1], 'box', 'uniform', 1e-3);
% >> exact = 0.25;
% >> success = double(abs(exact - integ) <= 1e-3)
% success = 1
%
%
% Example 2:
% Integral of exp(-x1^2-x1^2) over [-1 -1; 2 2] 
% >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); 
% >> absTol = 1e-3;
% >> relTol = 1e-2; 
% >> domain = [-1 -1; 2 2];
% >> integ = cubMC_CLT('f', f, 'absTol', absTol, 'relTol', relTol, 'domain', domain);
% >> exact = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
% >> success = double(abs(exact - integ) <= max(absTol,relTol*exact))
% success = 1



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
% 
% % % Example 6: 
% >> w.f = @(x) (x.^(-1/3.).^2);
% >> w.domain=[zeros(1,1); ones(1,1)];
% >> w.absTol = 1e-2; w.relTol=0;
% >> [q, out_param] = cubMC_CLT(w);
% exactsol=3;
% 
% % % Example 7:
% >> w.f= @(x) 1-(2*abs(x-0.5));
% >> w.domain=[zeros(1,1); ones(1,1)];
% >> w.absTol = 1e-3; w.relTol=0;
% >> [q, out_param] = cubMC_CLT(w);
% exactsol = 0.5;
%
% _Author_

t_start = tic;
cub_inp = gail.cubMCParam(varargin{:}); %parse the input and check it for errors
out = gail.cubMCOut(cub_inp); %create the output class
[integ, out] = meanMC_CLT(out);
out.time = toc(t_start);

end 