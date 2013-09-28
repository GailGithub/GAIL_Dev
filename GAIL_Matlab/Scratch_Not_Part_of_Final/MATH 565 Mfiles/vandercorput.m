%van der Corput sequence
function x=vandercorput(n)
[mantissa,m]=log2(n);
m=max(m,0);
twom=2^m;
xtemp=zeros(twom,1);
twoip1=1; twominip1=1;
for i=1:m;
    twoi=twoip1;
    twoip1=2*twoip1;
    twominip1=twominip1/2;
    %keyboard
    xtemp(twoi+1:twoip1)=xtemp(1:twoi)+twominip1;
end
x=xtemp(1:n);
%keyboard