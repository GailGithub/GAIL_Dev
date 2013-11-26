function [stded_hat,gamma_hat,et_hat]=EstimateStd(varargin)
% EstimateStd finds the estimation of the standard deviation of IPaS
% output. It also return other variables such as gamma_hat, the estimation of the
% true probability and time, the average time to calculate one estimation.
% Input parameters are the same one used as an input in the IPaS function.
%
%   stded_hat=EstimateStd(f,split,T,M) only returns the estimate of the
%   standard deviation for that particular inputs
%
%   [stded_hat,gamma_hat]=EstimateStd(coeff) returns the estimation of the standard
%   deviation of the estimation, and the estimation of the true value 
%
%   [stded_hat,gamma_hat,time]=EstimateStd(coeff) returns the estimation of the 
%   standard deviation of the estimation, the estimation of the true value
%   and the average time to estimate the true value using the given input
%   coeff
%
%   gamma_hat --- the estimation of gamma, the probability that the last
%   term of a sequence with T terms is greater or equal to a predefined
%   value "split(end)".
%
%   et --- estimation of elapsted time for the estimator ruled by coeff
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
% Example 1: Estimating standard deviation
% >> coeff.M=1000; coeff.split= [2 4 5];
% >> stded_hat = EstimateStd(coeff);
%
%
% Example 2: Estimating standard deviation and expected value and
% elapsed_time
% >> coeff.M=1000; coeff.split= [2 4 5];
% >> [stded_hat,gamma_hat,et_hat] = EstimateStd(coeff);
% >> gamma_hat
%
%gamma_hat =  0.001*** 
%
%
% See also, IPaS.m, Mutation.m
%
% [1] Silva, Tiago. 2013. Rare event simulation for financial modeling with
% Interacting Path Systems and quasi-Monte Carlo methods. MS Thesis.
% Illinois Institute of Technology.



parameters = IPaS_param(varargin{:});
x=zeros(100,2);
for i=1:100
[x(i,1),x(i,2)]=IPaS(parameters);
if(i==3); SendMessage(mean(x(1:3,2))); end;
end
stded_hat=std(x(:,1));
gamma_hat=mean(x(:,1));
et_hat=mean(x(:,2));

function SendMessage(tm) 
if(tm>0.9); disp(['The stdev estimation may take around  = ' num2str(tm*100/60,3) ' minutes...']); end;



