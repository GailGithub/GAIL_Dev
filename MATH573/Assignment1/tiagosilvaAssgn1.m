% Name: Tiago Silva
% Email: tsilva@hawk.iit.edu
%
% What I am expecting from this class is to learn how to proper
% write source codes for computer programs in a professional way, learn the 
% programming styles that are currently accpected and being used in the 
% market, learn how to create maintainable source codes that are
% accessible, justified, efficient and robust and learn the most major and
% minor errors/mistakes/badpractices that are commonly used by programmers,
% in order to avoid them in the future.

%% A successful call of funappx_g.m
f = @(x) x.*(1-x);
[fappx, out_param] = funappx_g(f,'ninit',10,'nmax',1e6,'abstol',1e-2)
x = 0:1/99:1;
plot(x,fappx(x))

%% Output of funappx_g.m
%
% fappx = 
%
%    @(x)interp1(x1,y1,x,'linear')
%
%
% out_param = 
%
%          abstol: 0.0100
%               f: @(x)x.*(1-x)
%           ninit: 10
%            nmax: 1000000
%             tau: 17
%    exceedbudget: 0
%         npoints: 37
%        errbound: 1.9290e-04
%

%% At least one item in the GAIL library that should be corrected or improved
% I believe that it could be interesting to include an IF statement at line
% 4 of GAIL_Install.m file checking whether the file GAILstart.m is at the
% same folder or not. Currently, if for any reason, the file GAILstart.m is 
% not there, a standard error will appers, like:
% Undefined function or variable 'GAILstart'.
% while it would be more interesting to have a specific obstervation 
% written by the modeler.
%
%clear all; close all; clc;
%if (exist('GAILstart','file')==2)
%  [GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
%  fprintf('\nWelcome to GAIL version %g.\n', GAILVERSION);
%  if MATLABVERSION < 7,
%    error('This version is only supported on Matlab 7.x and above.');
%  else
%    gailp=genpath(GAILPATH);% adding all subdirectories 
%  end
%  addpath(gailp);% Add GAIL directories and subdirectories
%  savepath; % Save the changes
%  fprintf('\nGAIL version %g has been installed successfully.\n\n', GAILVERSION);
%else
% fprintf('\nGAIL has not been installed. Please, verify if all files were unzipped properly \n\n');
%end
