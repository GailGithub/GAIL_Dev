function x=qeq(a,b,c)

% x=QUADRATICEQ(a,b,c) finds the roots of the quadratic equation
% a x^2 + b x + c = 0

x=NaN(2,1); %initialize roots

%scale the inputs
scale=max(abs([a b c]));
a1=a/scale;
b1=b/scale;
c1=c/scale;
if scale==0, return, end

%compute the roots
disc=sqrt(b1^2-4*a1*c1);
term=(-b1 - sign(b1)*disc);
if abs(term)>0
    x(1)=(2*c1)/(-b1 - sign(b1)*disc);
    x(2)=(-b1 - sign(b1)*disc)/(2*a1);
else
    x=zeros(2,1);
end
x=sort(x);
