function x=LanJiangqeq(varargin)
% x=QEQ(a,b,c) finds the roots of the quadratic equation
%   a x^2 + b x + c = 0

x=[]; %initialize roots
coeff = qeq_param(varargin{:});
%% scale the inputs
scale=max(abs([coeff.a coeff.b coeff.c]));
a1=coeff.a/scale; %scale coefficients to avoid overflow or underflow
b1=coeff.b/scale;
c1=coeff.c/scale;
if scale==0, return, end %zero polynomial
%% compute the roots
term=-b1 - sign(b1)*sqrt(b1^2-4*a1*c1); %no cancellation error here
if term~=0 % at least one root is nonzero
    x(1)=(2*c1)/term;
    if a1~=0; x(2)=term/(2*a1); end %second root exists
elseif a1~=0
    x=zeros(2,1);
end
x=sort(x);
end

function coeff = qeq_param(varargin)
default.a = 1;
default.b = 2;
default.c = 1;
p = inputParser;
if isempty(varargin)
    warning('The coefficient must be provided, now using default a b and c')
    coeff.a = default.a;
    coeff.b = default.b;
    coeff.c = default.c;
elseif (nargin == 1 && isstruct(varargin{1}))
    p.StructExpand = true;
    p.KeepUnmatched  = true;
    addParamValue(p,'a',default.a,@isnumeric);
    addParamValue(p,'b',default.b,@isnumeric);
    addParamValue(p,'c',default.c,@isnumeric);
    parse(p,varargin{:})
    coeff = p.Results;
elseif (nargin ==1 && isvector(varargin{1}) && length(varargin{1})==3)
    coeff.a = varargin{1}(1);
    coeff.b = varargin{1}(2);
    coeff.c = varargin{1}(3);   
elseif (nargin == 3 && isnumeric(varargin{1})) 
    addRequired(p,'a',@isnumeric);
    addRequired(p,'b',@isnumeric);
    addRequired(p,'c',@isnumeric);
    parse(p,varargin{:})
    coeff = p.Results;
else
    warning('your input could not be recognized, now using default setting.')
    coeff.a = default.a;
    coeff.b = default.b;
    coeff.c = default.c;
end

if coeff.a ==0 &&  coeff.b ~= 0
    warning(['a could not 0,since this is a quadratic equation, '...
        'now using default a'])
    coeff.a = default.a;
end
if  coeff.b == 0&& coeff.a == 0
        warning('a and b could not both be zero, now using default a and b')
        coeff.a = default.a;
        coeff.b = default.b;
end
end

% Example 1
% x = LanJiangqeq(1,2,1)
% 
% x =
% 
%     -1    -1


% Example 2
% x = LanJiangqeq([1,2,1])
% 
% x =
% 
%     -1    -1


% Example 3
% coeff.a = 1;
% coeff.b = 2;
% coeff.c = 1; 
% x = LanJiangqeq(coeff)
% 
% x =
% 
%     -1    -1


% Example 4
% x = LanJiangqeq(1)
%
% Warning: your input could not be recognized, now using
% default setting. 
% > In LanJiangqeq>qeq_param at 53
%   In LanJiangqeq at 6 
% 
% x =
% 
%     -1    -1



