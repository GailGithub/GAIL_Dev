function MVNPfunvalfinal = MVNPexact(t,b,sig)
% MVNPexact calculates the true solution of multivariate
% normal probability when the coveriance matrix is in a special form:
% diagnal is 1 and off diagnal are all same.
%
% b - the upper limits of the integal with size 1 x d
% sig - the off diagnal element
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