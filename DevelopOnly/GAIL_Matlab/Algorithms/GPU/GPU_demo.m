%% GPU Processing DEMO
% *Authors:* Fabio Araujo da Silva and Renan Guarese
%% Objective
% Improving the performance of GAIL routines through solutions in *Parallel
% Computing* and *GPU Processing*.
%% Methodology
% - Matlab's *Parallel Computing Toolbox (PCT)*: both CPU and GPU.
%% Materials and Methods
% - All tests were done using an 8-Core I7 2.6GHz CPU, with a NVIDIA GeForce GTX 965M GPU(1024 cores), on Windows 10; 
%%
% - Among the functions adapted, MeanMC_CLT was chosen to be displayed
% here;
%% 
% - *MeanMC_CLT:* Monte Carlo method toestimate the mean of a random
% variable;
%% Discussion 
% The _gpuArray_ library was used to increase the performance of the
% algorithms. Functions like _zeros_, _ones_ and _rand_ - including its
% variations _randn_ and _randi_ - can be generated in or coverted to
% _gpuArray_ structure.
%%
% Example 1:
%% 
% *This*
rand(3,1)
%%
% *Becomes*
gpuArray.rand(3,1)
%%
% Example 2:
%% 
% *This*
ones(3,1)
%%
% *Becomes*
gpuArray(ones(3,1))
%%
% The first example can be used on certain functions only, such as the ones
% cited above. The second example, on the other hand, can be used to
% convert any MATLAB array into gpuArray.
%%
% Once there is an gpuArray variable, all of the native functions (_var()_, _mean()_,_sum()_ ...),
% arithmetic and relational operators. 
%% Results: Numerical Examples
% All of the examples run meanMC_CLT.m (CPU) and meanMC_CLT_GPU.m (GPU), plotting a
% graph of the general performance of both algorithms. 
%% * Example 1
% $$ Max (S(T)-K,0),\quad where $$
%%
% $$ S(T) = S_0\times e^{(\frac{-\sigma^2}{2}\times T+\sigma\times
% \sqrt{\frac{T}{d}}\times \sum_{j=1}^{d} Z_j )}$$
%%
% $$ Z_j = IID\quad N(0,1) $$
X = [1,2,3];                         %executions of the program
Y = [1,2,3];                         %this vector will store the runtimes (CPU)
Y_GPU = [1,2,3];                     %this vector will store the runtimes (GPU)
S0 = 100;
sigma = 0.5;
d = 15;
T = 2;      
K = S0;
Yrand = @(n) max(S0*exp(((sigma*sigma)/2)*T+sigma*sqrt(T/d)*sum(randn(n,d),2)) - K ,0);    %new Yrand function (CPU)
Yrand_GPU = @(n) max(S0*exp(((sigma*sigma)/2)*T+sigma*sqrt(T/d)*sum(gpuArray.randn(n,d),2))  - K,0) ;   %new Yrand function (GPU)
for i = X
    %run on CPU
    t = tic;
    meanMC_CLT(Yrand,1e-3);          %function and error tolerance as parameters
    Y(i) = toc(t);
    
    %run on GPU
    t = tic;
    meanMC_CLT_GPU(Yrand_GPU,1e-3);      %function and error tolerance as parameters
    Y_GPU(i) = toc(t);
end
%plotting graph
plot(X,Y,X,Y_GPU);     
xlabel('executions');         
ylabel('time(s)');
title('Example 1');
legend('CPU','GPU');
set(gca,'xtick',1:3);
AverageSpeedup = mean(Y)/mean(Y_GPU)
%% * Example 2
% $$ Yrand = max(e^{Z_j},0),\quad where $$
%%
% $$ Z_j = IID\quad N(0,1) $$
X = [1,2,3];                         %executions of the program
Y = [1,2,3];                         %this vector will store the runtimes (CPU)
Y_GPU = [1,2,3];                     %this vector will store the runtimes (GPU)
Yrand = @(n) max(exp(randn(n,d)),0);        %new Yrand function (CPU)
Yrand_GPU = @(n) max(exp(gpuArray.randn(n,d)),0);     %new Yrand function (GPU)
for i = X
    %run on CPU
    t = tic;
    meanMC_CLT(Yrand,1e-3);          %function and error tolerance as parameters
    Y(i) = toc(t);
    
    %run on GPU
    t = tic;
    meanMC_CLT_GPU(Yrand_GPU,1e-3);      %function and error tolerance as parameters
    Y_GPU(i) = toc(t);

end
%plotting graph
plot(X,Y,X,Y_GPU);     
xlabel('executions');         
ylabel('time(s)');
title('Example 2');
legend('CPU','GPU');
set(gca,'xtick',1:3);
AverageSpeedup = mean(Y)/mean(Y_GPU)
%% * Example 3
% $$ Yrand = max(cot(\sum Z_j),0),\quad where $$
%%
% $$ Z_j = IID\quad N(0,1) $$
d = 3;
X = [1,2,3];                         %executions of the program
Y = [1,2,3];                         %this vector will store the runtimes (CPU)
Y_GPU = [1,2,3];                     %this vector will store the runtimes (GPU)
Yrand = @(n) max(cot(sum(rand(n,d),2)),0);        %new Yrand function (CPU)
Yrand_GPU = @(n) max(cot(sum(gpuArray.rand(n,d),2)),0);     %new Yrand function (GPU)
for i = X
    %run on CPU
    t = tic;
    meanMC_CLT(Yrand,1e-3);          %function and error tolerance as parameters
    Y(i) = toc(t);
    
    %run on GPU
    t = tic;
    meanMC_CLT_GPU(Yrand_GPU,1e-3);      %function and error tolerance as parameters
    Y_GPU(i) = toc(t);
end
%plotting graph
plot(X,Y,X,Y_GPU);     
xlabel('executions');         
ylabel('time(s)');
title('Example 3');
legend('CPU','GPU');
set(gca,'xtick',1:3);
AverageSpeedup = mean(Y)/mean(Y_GPU)
%% References
% [1] HICKERNELL, Fred J .et al. "GAIL: Guaranteed Automatic Integration Library", 2011. http://gailgithub.github.io/GAILDev/
%%
% [2] ALTMAN, Yair. "Accelerating MATLAB Performance", 2014: CRC Press
%%
% [3] REESE, Jill; ZARANEK, Sarah. "GPU Programming in MATLAB", 2012. http://www.mathworks.com/