function x=XuanZhouQeq(varargin)
% x=QEQ(a,b,c) finds the roots of the quadratic equation
%   a x^2 + b x + c = 0

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

function out_param=qeq_param(varargin)
%% Parse inputs
%p = inputParser;
if nargin == 3
    out_param.a = varargin{1};
    out_param.b = varargin{2};
    out_param.c = varargin{3};    
elseif isstruct(varargin{1})
    out_param = varargin{1};
elseif isvector(varargin{1})
    out_param.a = varargin{1}(1);
    out_param.b = varargin{1}(2);
    out_param.c = varargin{1}(3);
else
    error('The inputs should be of the form: XuanZhouQeq(a,b,c) or c]) or  XuanZhouQeq([a b c]) or  XuanZhouQeq(coef), where coef.a, coef.b, and coef.c are the coefficients.');
end

%% Check parameter validity
if not(isnumeric(out_param.a) && isnumeric(out_param.b) && isnumeric(out_param.c))
    error('The input coefficients should all be numeric.');
end
if out_param.a == 0
    error('The polynomial is linear or constant.');
end
end

% x = XuanZhouQeq(1,2,3)
% 
% x =
% 
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
% 
% x = XuanZhouQeq([1 2 3])
% 
% x =
% 
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
% 
% coef.a = 1; coef.b = 2; coef.c = 3;
% x = XuanZhouQeq(coef)
% 
% x =
% 
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
% 
% x = XuanZhouQeq(0,2,3)
% Error using XuanZhouQeq>qeq_param (line 47)
% The polynomial is linear or constant.
% 
% Error in XuanZhouQeq (line 6)
% out_param = qeq_param(varargin{:});