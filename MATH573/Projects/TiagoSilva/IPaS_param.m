function coeff_out = IPaS_param(varargin)
% IPaS_param validades the input and return a valid output than can be used
% as a input in another function. If input is not valid, it returns a
% valid default value.
%
%   coeff_out = IPaS_param(f,split,T,M)
%
%   coeff_out = IPaS_param(coeff_in) where coeff_in contains the following
%   parameters: coeff_in.f coeff_in.split, coeff_in.T, coeff_in.M
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
% Example1 : Empty Input
% >> IPaS_param()
%
% Warning: The coefficient must be provided, now using default input
% > In IPaS_param at 72 In doctest_run>DOCTEST__evalc at 57 In doctest_run at 25 In doctest at 174
% ans = 
%        f: @(v,U)v+(U<0.1)
%    split: [1 2 3 4 5 6]
%        T: 10
%        M: 1000
% 
%
% Example2 : Using only split levels
% >> coeff.split = [3,7];
% >> IPaS_param(coeff)
%
% ans = 
%        f: @(v,U)v+(U<0.1)
%        M: 1000
%    split: [3 7]
%        T: 10
%
%
% Example3 : Using ordered input.
% >> f = @(v,U) v+(U<0.1);split = [1,2,3,4,5,6]; T=10; M=10^4;
% >> IPaS_param(f,split,T,M)
%
% ans = 
%        f: @(v,U)v+(U<0.1)
%        M: 10000
%    split: [1 2 3 4 5 6]
%        T: 10
%
%
% See also, IPaS.m, Mutation.m
%
% Reference: MATH 573 Reliable Mathematical Software
%

%% Parsing
default.f = @(v,U) v+(U<0.1);
default.split = [1,2,3,4,5,6];
default.T = 10;
default.M=10^3;
p = inputParser;
if isempty(varargin)
    warning('The coefficient must be provided, now using default input')
    coeff.f = default.f;
    coeff.split = default.split;
    coeff.T = default.T;
    coeff.M = default.M;
elseif (nargin == 1 && isstruct(varargin{1}))
    p.StructExpand = true;
    p.KeepUnmatched  = true;
    addParamValue(p,'f',default.f,@isfcn);
    addParamValue(p,'split',default.split,@isvector);
    addParamValue(p,'T',default.T,@isnumeric);
    addParamValue(p,'M',default.M,@isnumeric);
    parse(p,varargin{:})
    coeff = p.Results;
 elseif (nargin == 4 && isfcn(varargin{1}) && isvector(varargin{2}))
    addRequired(p,'f',@isfcn);
    addRequired(p,'split',@isvector);
    addRequired(p,'T',@isnumeric);
    addRequired(p,'M',@isnumeric);
    parse(p,varargin{:})
    coeff = p.Results;
else
    warning('Your input could not be recognized, now using default setting.')
    coeff.f = default.f;
    coeff.split = default.split;
    coeff.T = default.T;
    coeff.M =default.M;
end

%% Validadtion
if ~isposint(coeff.T) % number of sample paths should be an integer
warning(['The number of steps should be a positive integer, using ' num2str(default.T)])
coeff.T = default.T;
end

if size(coeff.split,2)>1
    if(min(diff(coeff.split))<=0); % slipt levels must be an increasing sequence
    warning(['The levels defined on vector coeff.split must be a increasing sequence, '...
        'now using default [1,2,3,4,5,6]'])
    coeff.split = default.split;
    end
    
    if(size(coeff.split,2)>coeff.T) %number of split levels cannot be bigger than number of steps
    if (coeff.T>1)
       warning(['The number of levels must be less than the total number of steps coeff.T, '...
      'now using split = [1 2]'])
      coeff.split = [1 2];
    else
      warning(['The number of levels must be less than the total number of steps coeff.T, '...
      'now using split = [1]'])
      coeff.split = [1]; 
    end

    end
end

if (nargin(coeff.f)~=2) % function must accept two inputs (condition may change in the future)
     warning(['The function coeff.f must have two arguments, where the second one receives the random variable '...
       'now using default @(v,U) v+(U<0.1);'])
    coeff.f = default.f;
end
    
if ~isposint(coeff.M) % number of sample paths should be an integer
warning(['The number of samples should be a positive integer, using ' num2str(default.M)])
coeff.M = default.M;
end

coeff_out=coeff;

