function x=yizhiqeq(varargin)
% YIZHIQEQ finds the roots of the quadratic equation a x^2 + b x + c = 0.
% Yizhi Zhang yzhang97@hawk.iit.edu
% 
%  Description
%
%   x = YIZHIQEQ(a,b,c) computes x, roots of the quadratic equation a x^2 +
%   b x + c = 0. The input of a, b, and c has to be complete and in sequence. 
%
%   x = YIZHIQEQ(array) computes x, roots of the quadratic equation a x^2 +
%   b x + c = 0. The input is an array with three element in sequence of a, b, and c.
% 
%   x = YIZHIQEQ(coeff) computes x, roots of the quadratic equation a x^2 +
%   b x + c = 0. The input coeff is a structure that contains the parameter value a, b, and c in sequence. 
% 
%
%   Examples
%
%   Example 1: 
%   >> format short e; yizhiqeq(1,-5,4)
%    1.0000e+00   4.0000e+00
%
%
%   Example 2:
%   >> yizhiqeq([1,-5,4])
%   1.0000e+00   4.0000e+00
%
%
%   Example 3:
%   >>   coeff.a=1; coeff.b=-5; coeff.c=4; yizhiqeq(coeff)
%   1.0000e+00   4.0000e+00
%
%
%   Example 4:
%   >> q = yizhiqeq()
%   Warning: Coefficient must be specified. Now using default value a = 1, b = 3 and c = 2 
%   >  In ***
%   -2.0000e+00  -1.0000e+00
%
%
% See also yuhanqeq.m LanJiangqeq.m
%
% Reference:
% [1]  S. T. Choi, F. J. Hickernell, lecture slides MATH573slides.pdf

% check parameter satisfy conditions or not
coeff = yizhiqeq_param(varargin{:});

x=[]; %initialize roots
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

function coeff = yizhiqeq_param(varargin)
% parse the input to the integral_g function

% Default parameter values
default.a  = 1;
default.b  = 3;
default.c  = 2;


if isempty(varargin)
    warning(['Coefficient must be specified. Now using default value a = '...
        num2str(default.a) ', b = ' num2str(default.b) ' and c = ' num2str(default.c) ])
    coeff=default;
end;

validvarargin=numel(varargin)>=1;
if validvarargin
    in2=varargin{:};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if validvarargin
    p = inputParser;
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'c',default.c,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'c',default.c,@isnumeric);
    end
    parse(p,varargin{:})
    coeff = p.Results;
end;

if isvector(varargin) && length(varargin) == 1 && ~isstruct(in2)
    coeff.a = varargin{1}(1);
    coeff.b = varargin{1}(2);
    coeff.c = varargin{1}(3);
end;

if coeff.a ==0 &&  coeff.b ~= 0
    warning(' Polynomial is linear.')
end

if  coeff.b == 0&& coeff.a == 0
        warning('Polynomial is a constant.')
end

end


%% results
% doctests
% >> doctest yizhiqeq
% TAP version 13
% 1..4
% ok 1 - "format short e; yizhiqeq(1,-5,4)"
% ok 2 - "yizhiqeq([1,-5,4])"
% ok 3 - "  coeff.a=1; coeff.b=-5; coeff.c=4; yizhiqeq(coeff)"
% ok 4 - "q = yizhiqeq()"

