%This script will use a1=0.3 and a2=0.75 to construct a function with two
%  local minimum points.
%  
%  Generates Figure 3.2 in the thesis
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
a1=0.3;a2=0.75;
f=@(x) -5*exp(-(10*(x-a1)).^2)-exp(-(10*(x-a2)).^2);
t=0:0.00001:1;
plot(t,f(t),'LineWidth',2)
xlim([0 1])
ylim([-6 1])

save_eps('UniFunMinOutput', 'UniFunMinTwoExtreme');