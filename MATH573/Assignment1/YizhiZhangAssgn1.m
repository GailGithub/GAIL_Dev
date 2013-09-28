%Yizhi Zhang    yzhang97@hawk.iit.edu
%
%My goals for this course are 
%   o to learn more about writing numerical computation software in class
%   o to improve my programming skill by practising and discussing with
%   teachers and classmates
%   o to learn the implement of this software in finance problems.
% 

abstol=1e-8;
tic
z=0.5; a=0.1;
f=@(x) exp(-(c*(x-z)).^2).*(2*c./...
  (sqrt(pi)*(erf(c*(1-z))+erf(c*z)))); 

appxinteg=integral_g(f,'abstol',abstol,'nmax',1e8);
time=toc;
trueinteg=1;
error=trueinteg-appxinteg;
disp(['       true integral = ' num2str(trueinteg,14)])
disp(['approximate integral = ' num2str(appxinteg,14)])
disp(['               error = ' num2str(error,5)])
disp(['                 tol = ' num2str(abstol,5)])
disp(['          time taken = ' num2str(time,5) ' seconds'])

%  YizhiZhangAssgn1
%        true integral = 1
% approximate integral = 0.99999999993601
%                error = 6.3991e-11
%                  tol = 1e-08
%           time taken = 0.0053863 seconds

%Suggestion for improvement: since we talked about bring up the meat to the
%top of the algorithm, I suggest that we collapse the help documents by
%clicking the minus sign in front of the block. Along with moving the
%parsing to the end of the file, the meat will show up when people open the
%file.

