%This script will use a=0.03 and z=0.4 to construct an example of the bump
%  test functions.
% 
%  Generates Figure 3.1 in the thesis
%
%  Xin Tong, A Guaranteed, Adaptive, Automatic Algorithm for Univatiate
%  Function Minimization, July 2014.


%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
%Some defaults to make the plots easier to view
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultLineMarkerSize',40) %large dots

%% Plot
a = 0.03;z=0.4;
f=@(x) 0.5/a^2*(-4*a^2-(x-z).^2-(x-z-a).*abs(x-z-a)+(x-z+a).*...
    abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);
t=0:0.00001:1;
plot(t,f(t),'LineWidth',2)
xlim([0 1])
ylim([-1.2 0.2])
