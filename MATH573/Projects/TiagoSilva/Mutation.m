function [tau_o,y]=Mutation(f,x,tau_i,lvl,T,U)
% Mutation returns the value y, and the stopping time tau_o, given the 
% Markov transition function f:X_i->X_(i+1), the previous value x, the 
% previous stopping time tau_i, the goal level lvl, a finite time T and a 
% sequence of random variable U. Considering x as the (tau_i)^(th) term of 
% a particular sequence, Mutation apply the function f  j=1:(T-(tau_i)) 
% using the given random variable sequence U, until x >= lvl. In this 
% scenario, tau_o is the time j that the value x >= lvl and y=x at this 
% time j. If the value X_i does not reach the level lvl, i.e. X_T < lvl, 
% then we define the stopping time tau_o=(T+1).       
%
%   [tau_o,y]=Mutation(f,x,tau_i,lvl,T,U) finds the (tau_o)^(th) coefficient 
%   of the monotone Markov sequenc, where tau_o is the position where the 
%   coefficient is greater or equal the level lvl for the first time. All 
%   other inputs are described above.
%
%   [tau_o,y]=Mutation(f,x,tau_i,lvl,T,U) finds the (tau_o(i))^(th) coefficient 
%   of the monotone Markov sequenc, where tau_o(i) is the position where the 
%   coefficient is greater or equal the level lvl for the first time. All 
%   other inputs are described above. This function also accepts x and
%   tau_i as vectors. In this case, y and tau_o will also be vectors.
%
% Example 1: Using function to generate the parent for the first cycle.
% >> f = @(v,U) v+(U<0.1); T=10;U=[ ones(1,3) 0.05 ones(1,6)]; lvl=1;
% >> [tau_o,y]=Mutation(f,0,0,lvl,T,U)
%
% tau_o = 4
% y = 1 
%
%
% Example 2: Using function to generate the second stopping time from the
% previous cycle. Note that the vector U was set so that it reaches the
% second lvl=2 after three steps. Since the first stopping time is 4, then
% the next one will be 7.
% >> f = @(v,U) v+(U<0.1); T=10;U=[ones(1,2) 0.05 ones(1,6)];
% >> x=1; tau_i=4; lvl=2;
% >> tau_o=Mutation(f,x,tau_i,lvl,T,U)
%
% tau_o = 7
%
% 
% Example 3: Placing Inputs that is very unlikely that X_T >=lvl.
% >> f = @(v,U) v+(U<0.1); T=10;U=[ones(1,2)*0.04 0.3 ones(1,7)*0.07];
% >> x=0; tau_i=0; lvl=10;
% >> [tau_o,y]=Mutation(f,x,tau_i,lvl,T,U)
% tau_o = 11
% y = 9
% 
%
% See also, IPaS.m
%
% Reference: 
% [1] Del Moral, P., & Patras, F. (2011). Interacting path systems for 
% credit risk. In T. Bielecki, D. Brigo & F. Patras Credit Risk Frontier. 
% Bloomberg Press, Ch. 21, Sec. 4.
%
% [2] Silva, Tiago. 2013. Rare event simulation for financial modeling with
% Interacting Path Systems and quasi-Monte Carlo methods. MS Thesis.
% Illinois Institute of Technology.
%
N = size(x,1);          % Collect number of particles
tau=(T+1)*ones(N,1);    % Setting up the next stopping time
for i=1:N
    for j=1:(T-tau_i(i))
        % Generating next term of the monotone increasing sequence
        % using markov transition function f
        x(i)= f(x(i),U(i,j));   
        % If particle X reaches next level then collect stopping time
        if (x(i)>=lvl); tau(i) = j+tau_i(i);break; end;
    end
end
y=x;
tau_o=tau;

