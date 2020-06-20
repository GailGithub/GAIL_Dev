%% This is a plot of 256 random and Sobol points in [0,1)^2. 

%% Initialization
% initialize the workspace and set option parameters
clc;clearvars;

%% generate random and Sobol points
n = 2^8; %numbers of points to plot
d = 2;% dimension of points
rng(60616); %to get the same answer each time
sobolstr = sobolset(d);
xpts = [rand(n,d), sobolstr(1:n, 1:d)];% [256 randomly generated points, 256 Sobol points]

%% plot the points
for i = 1:2:3
figure
scatter(xpts(:,i),xpts(:,i+1), 25, 'k','filled')
switch(i)
   case 1
       title = './cvSobolPlot1.eps';
   case 3
       title = './cvSobolPlot2.eps';
end
print('-deps', title)
end
%
% %%
% _Author: Da Li 
