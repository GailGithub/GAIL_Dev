% Jagadeeswaran.R, jrathin1@iit.edu
%
% Windows8, Matlab 7.10.0 (R2010a)

function jagadeesHW3()
    % print help
    qeq
    
    % positive tests
    a=1; b=2; c=3;
	qeq(a,b,c)
	qeq([a b c])
    
    qeq('c', c, 'a', a, 'b', b)
    qeq('C', c, 'B', b, 'A', a)

	coef.a = a;
	coef.b = b;
	coef.c = c;
	qeq(coef)
    
    % warning
    a=0; b=2; c=3;
	qeq(a,b,c)
    a=0; b=0; c=3;
	qeq(a,b,c)

    a=1; b=0; c=3;
	qeq(a,b,c)

    a=Inf; b=0; c=3;
	qeq(a,b,c)

    % fail cases
    qeq(a,b)
	qeq([b c])
    qeq('c', c, 'a')
    qeq('X', c, 'Y', b, 'A', a)

	coef1.a = a;
	coef1.c = c;
	qeq(coef1)    
end

function x=qeq(varargin)
% x=QEQ(a,b,c) finds the roots of the quadratic equation
%   a x^2 + b x + c = 0

    % Parse parameters
	if isempty(varargin)
		fprintf('[x1 x2] = qeq(a, b, c) finds the roots of the quadratic equation\n');
		fprintf('   a x^2 + b x + c = 0\n');
		fprintf('[x1 x2] = qeq([a, b, c]) Input can be a vector\n');
		fprintf('[x1 x2] = qeq(''a'', a, ''b'', b, ''c'', c) Input can be lableled\n');
		fprintf('[x1 x2] = qeq(coeff) where coeff.a, coeff.b, coeff.c Input can be a structure\n');
		warning('Input arguments a,b,c must be specified.')
		return;
	end;
	
	n_argin = numel(varargin);
	
	if (n_argin>1)
		if isnumeric(varargin{1}) % format: qeq(a,b,c)
			if (n_argin==3) && isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
				a = varargin{1};
				b = varargin{2};
				c = varargin{3};
			else
				warning('Unknown Input arguments format.')
                return
			end
        elseif n_argin==6 % format ('b', b, 'a', a, 'c', c)
			
			for j=1:2:6
                switch varargin{j}
                    case {'a', 'A'}
                        a = varargin{j+1};
                    case {'b', 'B'}
                        b = varargin{j+1};
                    case {'c', 'C'}
                        c = varargin{j+1};
                    otherwise
                        warning('Unknown Input arguments format.')
                        return
                end
            end
        else
            warning('Unknown Input arguments format.')
            return
        end
	else
		if isstruct(varargin{1})
            v = varargin{1};
            if  isfield(v, 'a') && isfield(v, 'b') && isfield(v, 'c')
                a = v.a; b = v.b; c = v.c;
            else
                warning('Unknown Input arguments format.')
                return
            end
        elseif isvector(varargin{1}) && (numel(varargin{1})==3)
            v = varargin{1};
            a = v(1); b = v(2); c = v(3);
        else
            fprintf('Unknown Input arguments format.')
            return
		end
    end
    
    assert(exist('a','var') && exist('b','var') && exist('c','var') );

	x=[]; %initialize roots
	%% scale the inputs
	scale=max(abs([a b c]));
	a1=a/scale; %scale coefficients to avoid overflow or underflow
	b1=b/scale;
	c1=c/scale;
	if scale==0, return, end %zero polynomial
    if a1==0, warning('Linear polynomial.'), end
    if a1==0 && b1==0, warning('Constant polynomial.'), end
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

%>> jagadeesHW3
% [x1 x2] = qeq(a, b, c) finds the roots of the quadratic equation
%    a x^2 + b x + c = 0
% [x1 x2] = qeq([a, b, c]) Input can be a vector
% [x1 x2] = qeq('a', a, 'b', b, 'c', c) Input can be lableled
% [x1 x2] = qeq(coeff) where coeff.a, coeff.b, coeff.c Input can be a structure
% Warning: Input arguments a,b,c must be specified. 
%> In jagadeesHW3>qeq at 56
%   In jagadeesHW3 at 7
% 
% ans =
%
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
% 
% 
% ans =
%
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
% 
% 
% ans =
%
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
%
% 
% ans =
%
%   -1.0000 - 1.4142i  -1.0000 + 1.4142i
%
% 
% ans =
%
%  -1.0000 - 1.4142i  -1.0000 + 1.4142i
% 
%Warning: Linear polynomial. 
% > In jagadeesHW3>qeq at 118
%   In jagadeesHW3 at 24
%
% ans =
%
%    -1.5000
% 
%Warning: Linear polynomial. 
% > In jagadeesHW3>qeq at 118
%   In jagadeesHW3 at 26
%Warning: Constant polynomial. 
% > In jagadeesHW3>qeq at 119
%  In jagadeesHW3 at 26
% 
%ans =
% 
%      []
%
% 
% ans =
%
%      0
%      0
% 
%
% ans =
%
%    NaN   NaN
% 
% Warning: Unknown Input arguments format. 
% > In jagadeesHW3>qeq at 69
%   In jagadeesHW3 at 35
% Unknown Input arguments format.Warning: Unknown Input arguments format. 
% > In jagadeesHW3>qeq at 88
%   In jagadeesHW3 at 37
% Warning: Unknown Input arguments format. 
% > In jagadeesHW3>qeq at 83
%   In jagadeesHW3 at 38
%Warning: Unknown Input arguments format. 
% > In jagadeesHW3>qeq at 97
%  In jagadeesHW3 at 42
% >> 
% 

