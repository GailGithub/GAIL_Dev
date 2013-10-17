function x=silvaQEQ(varargin)
% Name: Tiago Silva
% Email: tsilva@hawk.iit.edu
%
%  The function silvaQEQ finds the roots of the quadratic equation
%  a x^2 + b x + c = 0
% 
%   The function accepts three different kinds of input:
%   silvaQEQ(a,b,c)
%   silvaQEQ([a,b,c])
%   silvaQEQ(coeff)  ; where coeff.a = a, coeff.b = b and coeff.c = c
%
%   and it returns a solution vector with the two roots
%   of the quadratic equation  a x^2 + b x + c = 0
%
%   EXAMPLES
%
%   Example 1
%   silvaQEQ([1,2,1])
%
%   ans =
%
%        -1    -1
%
%   Example 2
%   silvaQEQ(1,-2,1)
%   ans =
%
%         1     1
%
%   Example 3
%   coeff.a = 1; coeff.b = 3; coeff.c = 2; silvaQEQ(coeff)
%
%   ans =
%
%       -2.0000   -1.0000
%
%   See also: qeq.m
%
%   Reference: Course Reliable Mathematical Software
%              offered by Sou-Cheng T. Choi & Fred J. Hickernell

 


x=[]; %initialize roots
out_param = silva_param(varargin{:});

%% scale the inputs
scale=max(abs([out_param.a out_param.b out_param.c]));

if scale==0, return, end %zero polynomial

a1=out_param.a/scale; %scale coefficients to avoid overflow or underflow
b1=out_param.b/scale;
c1=out_param.c/scale;

if (a1==0&&b1~=0), x = -c1/b1; 
    warning(['Polinomial is linear']);
    return 
end
if (a1==0&&b1==0&&c1~=0)
    warning(['Polinomial is constant']);
    return 
end


%% compute the roots
term=-b1 - sign(b1)*sqrt(b1^2-4*a1*c1); %no cancellation error here
if term~=0 % at least one root is nonzero
    x(1)=(2*c1)/term;
    if a1~=0; x(2)=term/(2*a1); end %second root exists
elseif a1~=0
    x=zeros(1,2);
end
x=sort(x);
end


function out_param = silva_param(varargin)
default.a = 1;
default.b = -3;
default.c = 2;
valid = false;

% Collect nnumber of arguments
ni = nargin; 

if ni==3
    % If the case is silvaQEQ(a,b,c)
    if (isnumeric(varargin{1})&&isnumeric(varargin{2})&&isnumeric(varargin{3}))
        default.a = varargin{1};
        default.b = varargin{2};
        default.c = varargin{3};
        valid = true;
    end
    
elseif ni==1
 mainInput = varargin{1};

 if isstruct(mainInput)
    % If the case is silvaQEQ(coeff)
    default.a = mainInput.a;
    default.b = mainInput.b;
    default.c = mainInput.c;
    valid = true;

 elseif (isvector(mainInput)&&length(mainInput)==3)
    % If the case is silvaQEQ([a,b,c])
        default.a = mainInput(1);
        default.b = mainInput(2);
        default.c = mainInput(3);
        valid = true;
 end    
end


% Providing output
if (~valid)
    warning(['Input variable is not valid. Default variable will be used instead. a=1,b=-3,c=2']);
end
out_param.a = default.a;
out_param.b = default.b;
out_param.c = default.c;

end

%%OUTBPUTS

%Example 1
%silvaQEQ([1,2,1])
%
%ans =
%
%    -1    -1

%Example 2
%silvaQEQ(1,-2,1)
%ans =
%
%     1     1

%Example 3
%coeff.a=1; coeff.b=3; coeff.c=2; silvaQEQ(coeff)
%
%ans =
%
%   -2.0000   -1.0000

%Example 4
%silvaQEQ('a',1,'b',2,'c',2)
%Warning: Input variable is not valid. Default variable will be used instead.  a=1,b=-3,c=2  
%> In silvaQEQ>silva_param at 81
%  In silvaQEQ at 8 
%
%ans =
%
%    1.0000    2.0000