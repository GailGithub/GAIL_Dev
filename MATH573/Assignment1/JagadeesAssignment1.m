
% Jagadeeswaran R, jrathin1@iit.edu
%
% My goals for this course are
%   o To Improve my Programming skills for a Quality software development
%   o To understand techniques used in GAIL software development
%   o To contribute and learn from the team


%yf = @(x) exp(-abs(x));
%yf = @(x) x.*log(1+x.^2);
yf = @(x) (1/pi)*(1./(1+x.^2));
tic
[fappx, out_param] = funappx_g(yf);
time= toc;
x0 = 0.7658;
truefx0 = yf(x0);
fappxx0 = fappx(x0);
fprintf('True   f(%f) = %2.14f\n', x0, truefx0);
fprintf('Approx f(%f) = %2.14f\n', x0, fappxx0);
fprintf('Error        = %2.6e\n', abs(fappxx0-truefx0))
fprintf('Time taken   = %f seconds\n', time)


%  >> JagadeesAssignment1
% True   f(0.765800) = 0.20064291872750
% Approx f(0.765800) = 0.20064291953834
% Error        = 8.108454e-010
% Time taken   = 0.001667 seconds
 
% Suggestion: When i called fappx for x>1, I got NaN, Instead It may throw some warning
% 
