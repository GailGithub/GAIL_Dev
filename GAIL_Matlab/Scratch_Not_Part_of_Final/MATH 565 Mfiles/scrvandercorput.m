%scrambled van der Corput sequence
function x=scrvandercorput(n)
[mantissa,m]=log2(n);
m=max(m,0);
twom=2^m;
scrmat=eye(m); scrmat=scrmat+triu(rand(m,m)<0.5,1);
shvec=rand(1,m)<0.5;
xdig=zeros(twom,m);
twoip1=1;
for i=1:m;
    twoi=twoip1;
    twoip1=2*twoip1;
    xdig(twoi+1:twoip1,1:i)=[xdig(1:twoi,1:i-1) ones(twoi,1)];
end
xtemp=mod(xdig*scrmat+repmat(shvec,twom,1),2)*pow2((-1:-1:-m)');
%keyboard
xtemp=xtemp+rand(twom,1)/twoip1;
x=xtemp(1:n);
%keyboard