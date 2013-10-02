function x=qeq(a,b,c)
% x=QEQ(a,b,c) finds the roots of the quadratic equation
%   a x^2 + b x + c = 0

x=NaN(2,1); %initialize roots
%% scale the inputs
scale=max(abs([a b c]));
a1=a/scale; %scale coefficients to avoid overflow or underflow
b1=b/scale;
c1=c/scale;
if scale==0, return, end %zero polynomial
%% compute the roots
term=-b1 - sign(b1)*sqrt(b1^2-4*a1*c1); %no cancellation error here
if abs(term)>0 % at least one root is nonzero
    x(1)=(2*c1)/term;
    x(2)=term/(2*a1);
else
    x=zeros(2,1);
end
x=sort(x);
