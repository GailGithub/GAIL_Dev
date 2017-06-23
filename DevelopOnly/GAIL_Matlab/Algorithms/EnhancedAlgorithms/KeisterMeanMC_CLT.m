%% MeanMC_CLT & Keister's Example of Multidimensional Integration
% This MATLAB script shows how to use approximate Central Limit Theorem
% (CLT) with Monte Carlo to calculate Keister's Multidimensional Integration. 

%% Expressing the integral as an expectation
% Let's evaluate the integral using Monte Carlo cubature.  We first note
% that the change of variable \(\boldsymbol{t} = \boldsymbol{x}/a\)
% transforms this integral into
%
% \begin{align*} I &= \int_{\mathbb{R}^d} \cos(a \lVert \boldsymbol{t}
% \rVert) \exp(-a^2 \lVert \boldsymbol{t} \rVert^2) \, a^d \mathrm{d}
% \boldsymbol{t}, \qquad a > 0, \\ & = \int_{\mathbb{R}^d}
% \underbrace{(2\pi a^2)^{d/2} \cos(a \lVert \boldsymbol{t} \rVert)
% \exp((1/2-a^2) \lVert \boldsymbol{t} \rVert^2)}_{f(\boldsymbol{t})}
% \times \underbrace{\frac{\exp(-\lVert \boldsymbol{t} \rVert^2/2)}
% {(2\pi)^{d/2}}}_{\varrho(\boldsymbol{t})} \, \mathrm{d} \boldsymbol{t} \\
% & = \mathbb{E}[f(\boldsymbol{T})], \qquad \text{where } \boldsymbol{T} \sim \mathcal{N}(\boldsymbol{0},
% \mathsf{I}). \end{align*}

%% Evaluating the integral using |meanMC_g|
% To find \(I\) by Monte Carlo methods we define an anonymous function
% \(f\) as follows:

gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);

%%
% Next we call |meanMC_CLT|

abstol = 0; %absolute error tolerance
reltol = 0.01; %relative error tolerance
dvec = 1:5; %vector of dimensions
avec = [1 1/sqrt(2)]; %default value of a 
IMCvec = zeros(size(dvec)); %vector of answers
f2= @(t,d) cell2mat(arrayfun(@(a) f(t,a,d),avec,'UniformOutput',false));
tic
 for d = dvec
 f3=@(t)f2(t,d);
 YXn=@(n)f3(randn(n,d));
 s=struct('Y',YXn,'nY',size(avec,2)); %create a structure as input for meanMC_CLT
 IMCvec(d)= meanMC_CLTKATE(s,abstol,reltol);
 end

toc
%% Checking the real error
% There is a way to get the value of this integral to machine precision
% using the function |Keistertrue|
%
% <include>Keistertrue.m</include>

[~,Ivec] = Keistertrue(dvec(end));
relErrMC = abs(Ivec-IMCvec)./abs(Ivec);


