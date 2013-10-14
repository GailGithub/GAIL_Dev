function x=Xinqeq(varargin)
% x=Xinqeq(a,b,c) finds the roots of the quadratic equatio  a x^2 + b x + c = 0
%
% Xin Tong xtong5@hawk.iit.edu
%
%
%   x=Xinqeq(a,b,c) finds the roots of the quadratic equation a x^2 + b x + c = 0
%
%   x=Xinqeq([a,b,c]) finds the roots of the quadratic equation a x^2 + b x +
%   c = 0 whose coefficients are inputted by a vector.
%
%   x=Xinqeq(coef) finds the roots of the quadratic equation a x^2 + b x +
%   c = 0 whose coefficients are inputted by a structure coef with coef.a=a,
%   coef.b=b and coef.c=c.
%
%
%  Examples
%
%  Example 1:
%
%
%  >> x=Xinqeq([1,4,3])
% 
%   x =
% 
%     -3    -1
%
% 
% 
%  Example 2:
% 
% 
%  >>x=Xinqeq(1,3,2)
% 
%   x =
% 
%    -2.0000   -1.0000
% 
% 
% 
%  Example 3:
%  >>coef.a=1;coef.b=13;coef.c=4;
%  >>x=Xinqeq(coef)
% 
%  x =
% 
%   -12.6847   -0.3153
%
%
% 
%  See also
%
%  Reference



%% Parse and check the validity of input parameters
out_sample=Xinqeq_param(varargin{:});

x=[]; %initialize roots

%% scale the inputs
scale=max(abs([out_sample.a out_sample.b out_sample.c]));
a1=out_sample.a/scale; %scale coefficients to avoid overflow or underflow
b1=out_sample.b/scale;
c1=out_sample.c/scale;
if scale==0, return, end %zero polynomial

%% compute the roots
term=-b1 - sign(b1).*sqrt(b1.^2-4.*a1.*c1); %no cancellation error here
if term~=0 % at least one root is nonzero
    x(1)=(2*c1)./term;
    if a1~=0; x(2)=term./(2*a1); end %second root exists
elseif a1~=0
    x=zeros(2,1);
end
x=sort(x);
end



function out_sample=Xinqeq_param(varargin) % parse the input to the Xinqeq.m function
default.a=1;
default.b=3;
default.c=2;

%% Parse inputs
varnum=nargin;
in1=varargin{1};
if varnum==1 && isvector(in1) && length(in1)==3
    out_sample.a=in1(1);
    out_sample.b=in1(2);
    out_sample.c=in1(3);
else
    if  (varnum==1 && isstruct(in1))  || varnum==3
        p = inputParser;
        if isnumeric(in1)
            addOptional(p,'a',default.a,@isnumeric);
            addOptional(p,'b',default.b,@isnumeric);
            addOptional(p,'c',default.c,@isnumeric);
        else
            if isstruct(in1)
                p.StructExpand = true;
                p.KeepUnmatched = true;
            end
            addParamValue(p,'a',default.a,@isnumeric);
            addParamValue(p,'b',default.b,@isnumeric);
            addParamValue(p,'c',default.c,@isnumeric);
        end
        parse(p,varargin{1:end});
        out_sample = p.Results;
    else
        help Xinqeq
        warning(['You must input the three coefficients with the default three forms. Now Xinqeq is using a = ' num2str(default.a) '  b = ' num2str(default.b)  '  c = ' num2str(default.c)]);
        out_sample.a=default.a;
        out_sample.b=default.b;
        out_sample.c=default.c;
    end
end
 
%% Check parameter validity
if out_sample.a==0
    if out_sample.b==0
        warning(['The polynomial is constant.']);
    else 
        warning(['The polynomial is linear.']);
    end 
    out_sample.a=default.a;
    out_sample.b=default.b;
    out_sample.c=default.c;
end

end 






