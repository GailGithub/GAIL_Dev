%Yuhan Ding, yding2@hawk.iit.edu
%
%My goals for this course are 
%   o to increase the reliability and efficency of funappx in GAIL
%   o to gain more ideas about how to develop a reliable software
%   o to create more user-friendly documentation of GAIL

f = @(x) exp(sin(x));
tic
[fappx, out_param] = funappx_g(f);
time= toc;
x0 = sqrt(2)/2;
truefx0 = f(x0);
fappxx0 = fappx(x0);
disp(['true f(' num2str(x0) ') = ' num2str(truefx0,14)])
disp(['approximate f(' num2str(x0) ') = ' num2str(fappxx0,14)])
disp(['error = ' num2str(abs(fappxx0-truefx0),6)])
disp(['time taken = ' num2str(time,6) ' seconds'])

% YuhanDingAssgn1
% true f(0.70711) = 1.9148454972299
% approximate f(0.70711) = 1.9148454967427
% error = 4.87139e-10
% time taken = 0.00360072 seconds

%Suggestion for improvement: For unit/doc test, maybe a format or a
%consistent pattern for different algorithms is better for new developer to
%understand and implement.
