function MVNPfunvalfinal = MVNPexact(t,b,sig)
% this is the integrand used to calculate the true solution of multivariate
% normal probability when the coveriance matrix is a in a special form:
% diagnal is 1 and off diagnal are all the same.
%
% b is the upper limits of the integal with size 1 x d
% sig is the off diagnal element
% dim is the dimension of the integral
% t is the variable
MVNPfunval = (gail.stdnormcdf((b(1)+sqrt(sig)*t)/sqrt(1-sig)));
dim =  length(b);
for i =2:dim
    MVNPfunval= MVNPfunval.*(gail.stdnormcdf((b(i)+sqrt(sig)*t)/sqrt(1-sig)));
    %i=i+100;
end

MVNPfunvalfinal = MVNPfunval.*exp(-t.^2/2);

end