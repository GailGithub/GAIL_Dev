%% GUI of Funappx_g 
% To approximate a peaky function to show how *funappx_g* generate number
% of points

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
% We use *funappx_g* to approximate \(f\) over the interval \([0,1]\) with
% error tolerance \(1e-2\) and 15 initial subintervals
[~,out_param] = funappx_g(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,15)

% We can find to reach the error tolerance, we need 105 points to
% approximate the function. 

%% Process to Generate Number of Points
%
% Step 1: 16 evenly spaced points
%
% <<localgui1.png>>
% 
% Step 2: add points at peaky part
%
% <<localgui2.png>>
% 
% Step 6: after serveral iterations
%
% <<localgui6.png>>
% 
% Step 7: reach the error tolerance
%
% <<localgui7.png>>
% 
%
% This process can be reproduced by following command:
% funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,15,15);
%
%% References
%  
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%     _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%     Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%     GAIL: Guaranteed Automatic Integration Library (Version 2.2) [MATLAB
%     Software], 2017. Available from <http://gailgithub.github.io/GAIL_Dev/
%     GitHub>.

