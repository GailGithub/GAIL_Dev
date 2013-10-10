function x=yuhanqeq(varargin)
%YUHANQEQ Finds the roots of the quadratic equation a x^2 + b x + c = 0
%
%   x = YUHANQEQ(a,b,c) for given the ordered input parameters a,b,c,
%   ruturns the roots of a x^2 + b x + c = 0
%
%   x = YUHANQEQ([a,b,c]) or x = YUHANQEQ([a;b;c]) for given the ordered
%   input parameter a,b,c as a row vector or column vector, returns the
%   roots of a x^2 + b x + c = 0
%   
%   x = YUHANQEQ(coef) returns roots of the quadratic equation ax^2 + bx +c
%   = 0 where a = coef.a, b = coef.b and c = coef.c. If a field is not
%   specified, the default value is used.
%
%   x = YUHANQEQ('a',a,'b',c,'c')returns roots of the quadratic equation
%   ax^2 + bx +c = 0. All three field-value pairs are optional and can be
%   supplied in different order.
%
%   Examples
%
%   Example 1:
%   
%   >> x = yuhanqeq
%   
%   Warning: The coefficients of the quadratic equation must be
%   provided.Now using default a = 1 , b = -1 , and c = -2 > In
%   yuhanqeq>qeq_param at 138 In yuhanqeq at 100
%
%   x =
%   -1     2       
%
%
%   Example 2:
%   >> x = yuhanqeq(2,3,1)
%   
%   x =
%  -1.0000e+00  -5.0000e-01
%
%
%   Example 3:
%   
%   >> x = yuhanqeq([2,3,1])
%   
%   x =
%  -1.0000e+00  -5.0000e-01
%
%
%   Example 4:
%   
%   >> x = yuhanqeq([2;3;1])
%   
%   x =
%  -1.0000e+00  -5.0000e-01
%
%
%   Example 5:
%   
%   >> coef.a = 2; coef.b = 3; coef.c =1; x = yuhanqeq(coef)
%   
%   x =
%  -1.0000e+00  -5.0000e-01
%
%
%   Example 6:
%   
%   >> x = yuhanqeq('a',2,'b',3,'c',1)
%   
%   x =
%  -1.0000e+00  -5.0000e-01
%
%
%   Example 7:
%
%   >> x = yuhanqeq(0,1,3)
%
%   Warning:  coefficient of quadratic term is zero, quadratic equation
%   becomes a linear equation. > In yuhanqeq>qeq_param at 173 In yuhanqeq
%   at 100
% 
%   x =
%        -3
%
%
%   Example 8:
%
%   >> x = yuhanqeq(0,0,3)
%
%   Warning: a and b could not be zero, no root can be found > In
%   yuhanqeq>qeq_param at 179
%    In yuhanqeq at 100
% 
%    x =
% 
%          []
%
%
% YuhanDing yding2@hawk.iit.edu

x=[]; %initialize roots
out_param = qeq_param(varargin{:});

%% scale the inputs
scale=max(abs([out_param.a out_param.b out_param.c]));
a1=out_param.a/scale; %scale coefficients to avoid overflow or underflow
b1=out_param.b/scale;
c1=out_param.c/scale;
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

function out_param = qeq_param(varargin)
% parse the input to the yuhanqeq function

%% Default parameter values

default.a = 1;
default.b = -1;
default.c = -2;


validvarargin=numel(varargin)>=1;
if validvarargin
    in2=varargin{:};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

%if there is no input, choose default coefficents
if ~validvarargin
    warning(['The coefficients of the quadratic equation must be provided.'...
        'Now using default a = ' num2str(default.a) ' , b = ' num2str(default.b)...
        ' , and c = ' num2str(default.c) ])
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.c = default.c;
else
    p = inputParser;
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'c',default.c,@isnumeric);
    else
        if isstruct(in2)  %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end;
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'c',default.c,@isnumeric);
    end
    parse(p,varargin{:})
    out_param  = p.Results;
end

%if input is a vector
if isvector(varargin) && length(varargin) == 1 && ~isstruct(in2)
    out_param.a = varargin{1}(1);
    out_param.b = varargin{1}(2);
    out_param.c = varargin{1}(3);
end;

%if coefficient of quadratic term is zero
if out_param.a ==0 &&  out_param.b ~= 0
    warning([' coefficient of quadratic term is zero, quadratic equation '...
        'becomes a linear equation.'])
end

%if coefficients of quadratic term and linear term are zero
if  out_param.b == 0&& out_param.a == 0
        warning('a and b could not be zero, no root can be found')
end
end





