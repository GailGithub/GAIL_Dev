%Xuan Zhou, xzhou23@hawk.iit.edu
%
%My goals for this course are 
%   o to learn how to use Matlab to develop a numerical software of a
%   medium or large scale
%   o to gain some insight in the GAIL  
%   o to show support this new course as well as its instructors.

f = @(x) sin(x);
tic
[fappx, out_param] = funappx_g(f);
time=toc;
x0 = .5;
truefx0 = f(x0);
fappxx0 = fappx(x0);
disp(['true f(' num2str(x0) ') = ' num2str(truefx0,14)])
disp(['approximate f(' num2str(x0) ') = ' num2str(fappxx0,14)])
disp(['error = ' num2str(fappxx0-truefx0,5)])
disp(['time taken = ' num2str(time,5) ' seconds'])

% XuanZhouAssgn1
%true f(0.5) = 0.4794255386042
%approximate f(0.5) = 0.4794255386042
%error = 0
%time taken = 0.0033975 seconds

%Suggestion for improvement: One line in the documentation of MeanMC_g reads:
%in_param.npcmax --- number of elements in an array of optimal size to calculate the mu, the default value is 1e6.
%I recommended to change "the mu" to "mu" or "the mean".
