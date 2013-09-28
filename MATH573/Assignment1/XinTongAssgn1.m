% Xin Tong, xtong5@hawk.iit.edu
%
%My goals for this courses are
%   1. to learn what is a good software and how to write perfect codes
%   2. to meet new friends to learn more and get different kinds of information 

abstol=1e-12;
tic
f = @(x) -exp(-10000*(x-0.5).^2);
[fappx, out_param] = funappx_g(f,'ninit',10,'abstol',1e-8);
time=toc;
disp(['       tolerance = ' num2str(out_param.abstol)]);
disp(['        funcappx = ' char(@(x) interp1(x1,y1,x,'linear'))]);
disp(['         npoints = ' num2str(out_param.npoints)]);
disp(['        errbound = ' num2str(out_param.errbound)]);
disp(['      time taken = ' num2str(time,5) ' seconds']);


% Xin Tong  Assgn1
%         tolerance = 1e-08
%         funcappx = @(x)interp1(x1,y1,x,'linear')
%            npoints = 10
%         errbound = 4.929e-15
%       time taken = 0.0066324 seconds

%    This is a special peaky function f(x)= @(x) -exp(-10000*(x-0.5)^2) with funappx_g. 
%    If I choose a verry small initial number of points (ninit=10),  then the result shows 
%    that the we can just use 10 points to approximate the original function with our tolerance 1e-8 without
%    the warning "Warning: This function is peaky relative to ninit. You may wish to increase ninit
%    for similiar functions. " 

% Suggestion for improvement: to find an available initial n at first. 


