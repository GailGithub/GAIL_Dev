%% A GUI (graphical user interface) for *funappx_g*
% Author: Yuhan Ding, July 2017
%
% To approximate a peaky function with *funappx_g* and to show how
% *funappx_g* generates grid points for locally adaptive linear spline
% approximation

%% Function definition
%
% Define a peaky function as follows:
%
close all; clear all; format compact; format short;
f = @(x) exp(-1000*(x-0.2).^2);
x = 0:0.0001:1;
figure;
plot(x,f(x))
axis tight

%% Function Approximation
% We use *funappx_g* to approximate $f$ over the interval $[0,1]$ with
% error tolerance $10^{-2}$ and 15 initial subintervals:
[~,out_param] = funappx_g(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,15)

% We find that to reach the error tolerance, we need 105 points to
% approximate the function.

%% Process to Generate Grid Points
%
% Step 1: start with $16$ evenly spaced points:
%
% <<localgui1.png>>
%
% Step 2: add points to the peaky part:
%
% <<localgui2.png>>
%
% Step 6: after several iterations, the approximation error almost meets the given tolerance:
%
% <<localgui6.png>>
%
% Step 7: the error tolerance is reached:
%
% <<localgui7.png>>
%
%
% This process can also be reproduced by the following command:
% funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,15,15);
%
%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%     _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
