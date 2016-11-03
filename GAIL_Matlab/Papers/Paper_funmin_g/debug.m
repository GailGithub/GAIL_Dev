m=20;
s=zeros(1,m)
for i=1:m
    [fmin, fmax, out_min, out_max]=fejer_jackson_inequality(5*i);
    s(i)=size(out_min.intervals,2);
end
s

   