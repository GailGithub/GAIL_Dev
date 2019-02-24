function MVNPfunvalfinal = MVNPexact(t,b,sig)
% MVNPexact calculates the true solution of multivariate
% normal probability when the covariance matrix is in a special form:
% diagonal is 1 and off diagonal are all same.
%
% b - the upper limits of the integral with size 1 x d
% sig - the off diagonal element
% dim- the dimension of the integral
% t - the variable
MVNPfunval = (gail.stdnormcdf((b(1)+sqrt(sig)*t)/sqrt(1-sig)));
dim =  length(b);
for i =2:dim
    MVNPfunval= MVNPfunval.*(gail.stdnormcdf((b(i)+sqrt(sig)*t)/sqrt(1-sig)));
    %i=i+100;
end

MVNPfunvalfinal = MVNPfunval.*exp(-t.^2/2);

end