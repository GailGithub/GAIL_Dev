%inverse Gamma function
function y=invgam(x,a,b)
if nargin < 3; b=1;
    if nargin < 2; a = 1; end
end
myfun=@(y)gammainc(exp(y),a);
y=x;
x=x(:); %make x into a vector
[xsort,wh]=sort(x); %sort the x to make lookup easier
n=length(x);

%compute some values and perform linear interpolation
ysmall=0; xsmall=myfun(ysmall);
while xsmall>0.01; ysmall=ysmall-1; xsmall=myfun(ysmall); end 
ybig=0; xbig=myfun(ybig);
while xbig<0.99; ybig=ybig+1; xbig=myfun(ybig); end 
ytrial=ysmall+(ybig-ysmall)*(0:0.01:1);
xtrial=myfun(ytrial);
y0=interp1(xtrial,ytrial,xsort,'linear','extrap'); %interpolation for first guess

%perform secant method
y1=[y0(2:n); 2*y0(n)-y0(n-1)];
f0=myfun(y0);
f1=[f0(2:n); myfun(y1(n))];
f0=f0-xsort;
f1=f1-xsort;
ok=0;
%keyboard
while ok==0;
    y2=y1;
    whok=(f1~=f0); %protect against round-off error
    y2(whok)=y1(whok)-f1(whok).*(y1(whok)-y0(whok))./(f1(whok)-f0(whok));
    y0=y1; f0=f1;
    y1=y2; f1=myfun(y1)-xsort;
    err=max(abs(f1));
    if err<1e-12; ok=1; end
    %keyboard
end
y(wh)=b*exp(y0);