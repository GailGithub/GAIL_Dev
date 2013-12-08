function [gamma]=FindExactSolForBinoProblem(k,T,p)
% FindExactSolForBinoProblem finds the exact value for the probability that
% the sum of success in T bernoulli trials, which probability of success
% equal p, is greater or equal than k.
%
% [gamma] = FindExactSolForBinoProblem(k,T,p) returns the exact value
% gamma, where gamma = P(n#success >= k)
%
% gamma --- exact value for the given problem
%
% k --- minimum number of acceptable number of success. 
%
% T --- Total number of trials
%
% p --- probability of success for each trial
%
% Example 1: Simple test
% >> gamma = FindExactSolForBinoProblem(1,10,0.1)
% gamma = 0.6513
%
%
% Example 2: Extreme test
% >> gamma = FindExactSolForBinoProblem(7,10,0.1)
% gamma = 9.1216e-06
%
%
% See also, workout_1.m, workout_2.m

value=0;
for i=k:T
value = value + nchoosek(T,i)*p^i*(1-p)^(T-i);
end
gamma=value;