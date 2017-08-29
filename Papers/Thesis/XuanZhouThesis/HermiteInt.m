function result = HermiteInt(J)
H = HermitePoly(J);
result = 0;
for j = 1:J+1
    result = result+2*H(j)/j;
end
