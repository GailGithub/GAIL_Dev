%% Guarantee
%    
% If the function to be integrated, \(f\) satisfies the cone condition
%
% \[\|f''\|_1 \le \frac { \mathrm{nstar} }{2(b-a)}
% \left\|f'-\frac{f(b)-f(a)}{b-a}\right\|_1,\]
% 
% then the \(q\) output by this algorithm is guaranteed to satisfy
%
% \[\left\| \int_{a}^{b} f(x) dx - q \right\|_{1} \le \mathrm{abstol},\]
%
% provided the flag \(\mathrm{exceedbudget} = 0.\)
%
% And the upper bound of the cost is
%
% \[\sqrt{ \frac{\mathrm{nstar}* (b-a)^2 \mathrm{Var}(f')}{2 \times \mathrm{abstol}}}
% + 2 \times \mathrm{nstar} +4.\]
%
%
%% Examples
% *Example 1*

f = @(x) x.^2; [q, out_param] = integral_g(f)

% Integrate function x with default input parameter to make the error less
% than 1e-7.

%%
% *Example 2*

[q, out_param] = integral_g(@(x) exp(-x.^2),'a',1,'b',2,...
   'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)

% Integrate function x^2 with starting number of points 52, cost budget
% 10000000 and error tolerance 1e-8
