function [gamma,elapsed_time]=IPaS(varargin)
% IPaS finds the estimation of probability that X_T is greater or equal to 
% split(end), where the terms of the sequence X=(X_i), i=1..T, is ruled by 
% a Markov transition represented by the given function f:X_i -> X_(i+1), 
% and the vector split contains the levels for each cycle of the IPaS 
% method. The quantity M is the number of sample size. It is considered 
% that each term of X is a uni-dimensional random variable and the sequence 
% X is consequently T-dimentsional random variable.
%
%   [gamma_hat,et]=IPaS(f,split,T,M) finds the estimation of gamma, which 
%   is the probability that the last term X_T of a non decreasing sequence 
%   X ruled by the markov transition function f reaches the value 
%   split(end) using a sample size M through IPaS technique. Function also 
%   returns the elapsed_time et.
%
%   [gamma_hat,et]=IPaS(coeff) finds the estimation of gamma, which is the 
%   probability that the last term X_(coeff.T) of a non decreasing sequence 
%   X ruled by the markov transition function coeff.f reaches the value 
%   coeff.split(end) using a sample size coeff.M through IPaS technique. 
%   Function also returns the elapsed_time et.
%
%   gamma --- true value that we are trying to estimate.
%
%   gamma_hat --- the estimation of gamma, the probability that the last
%   term of a sequence with T terms is greater or equal to a predefined
%   value "split(end)".
%
%   et --- elapsted time for the estimator ruled by coeff
%
%   coeff.f --- a nondecreasing function f:(x,U)->y that is used to
%   generate the markov sequence X=(X_1,X_2,...X_T). the function f must
%   have two arguments, where the second one receives the random variable, 
%   the first one receives one term of the markvov sequence and the output
%   generate the next term of the same sequence.
%
%   coeff.split --- a vector that contais all the levels used for the IPaS
%   technique. If the vector contains only one coefficient, then the IPaS
%   method converge to naive MC method.
%
%   coeff.T --- total number of time steps and also the size of the markov
%   sequence X
% 
%   coeff.M --- sample size used to estimate gamma
%
% Example 1: Using ordered input.
% >> f = @(v,U) v+(U<0.1); split = [2,4,5]; T=10; M=10^4;
% >> gamma=IPaS(f,split,T,M)
%
% gamma = 0.001***
%
%
% Example 2: Using the structured input.
% >> coeff.f = @(v,U) v+(U<0.1); coeff.split = [2,4,5]; 
% >> coeff.T=10; coeff.M=10^4;
% >> gamma = IPaS(coeff)
%
% gamma = 0.001***
%
%
% See also, Mutation.m
%
% Reference: 
% [1] Del Moral, P., & Patras, F. (2011). Interacting path systems for 
% credit risk. In T. Bielecki, D. Brigo & F. Patras Credit Risk Frontier. 
% Bloomberg Press, Ch. 21, Sec. 4.
%
% [2] MATH 573 Reliable Mathematical Software
%
% [3] Silva, Tiago. 2013. Rare event simulation for financial modeling with
% Interacting Path Systems and quasi-Monte Carlo methods. MS Thesis.
% Illinois Institute of Technology.
%
tic;
output_param = IPaS_param(varargin{:});
f=output_param.f;
split=output_param.split;
T=output_param.T;
M=output_param.M;


x=zeros(M,1);           % Defining initial position for particles
t=zeros(M,1);           % Defining initial stopping time
kappa=size(split,2);    % Number of cycles for IPaS
m=ones(kappa,1);        % m(i) -- survival particles for each cycle
X=rand(M,T);            % Generating random variable for first cycle

[t,x]=Mutation(f,x,t,split(1),T,X); % S_1 - first set of particles xi
u=(t<=T);                           % xi that reached the first lvl
v=cumsum(u,1);                      % cummulative sum of vector u
m(1)= v(end);                       % Number of xi that reached the first lvl

%% Now, for each cycle we have first a Selection Phase and a Mutation Phase
for k=2:kappa
  % If zero particles reached the next level then there is nothing to be 
  % selected and IPaS is truncated here and returns gamma_hat=0  
  if (v(end)==0); m(k:kappa)=0; break; end;  
%% Selection Phase - Selecting only the particles that reached the 
%  previous level and replacing the others randomly selecting particles that
%  reached the previous level
  for i=1:M
    if(~u(i))
      j=find(v>rand*m(k-1),1); % randomly find a selectable index
      x(i)=x(j); % replace particle position
      t(i)=t(j); % replace last stopping time    
    end
  end
 % At this point, we have a new set where all particles 
 % reaches the previous level. This new set is called \tilda(S_(k-1))
 
%% Mutation Phase - Using this new set, \tilda(S_(k-1)) we are going to 
% generate a new set, S_(k) accoring to the given function.
X=rand(M,T-min(t));                 % Generate new random variable
[t,x]=Mutation(f,x,t,split(k),T,X); % Perform Mutation operator
u=(t<=T);                           % xi that reached the next lvl
v=cumsum(u,1);                      % cummulative sum of vector u
m(k)= v(end);                       % Number of xi that reached next lvl
end

gamma=prod(m/M);                    % estimating gamma from cond. prob.
elapsed_time=toc;                   % return estimation elapsed_time 



    
    
    
