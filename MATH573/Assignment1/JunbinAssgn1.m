%Juninbin Li, junbin731@gmail.com
%
%My goals for this assignment is:
%   o to know how to install run the MATLAB toolbox, Gail_Dev
%   o recognize the function in the software
%   

f = @(x) sqrt(cos(x)+1)  
tic
[fappx, out_param] = funappx_g(f);
time = toc;
x0 = 0.8;
truefx0 = f(x0);
fappxx0 = fappx(x0);
disp(['true f(' num2str(x0) ') = ' num2str(truefx0,14)])
disp(['approximate f(' num2str(x0) ') = ' num2str(fappxx0,14)])
disp(['error = ' num2str(fappxx0-truefx0,5)])
disp(['time taken = ' num2str(time,5) ' seconds'])


% true f(0.8) = 1.3025769494917
% approximate f(0.8) = 1.3025769450256
% error = -4.4661e-009
% time taken = 0.0014667 seconds

% I will add this part next time.

