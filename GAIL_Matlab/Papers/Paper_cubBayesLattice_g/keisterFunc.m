
function y = keisterFunc(x,dim,a)
  % a = 0.8;
  normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
  yinv = @(t) gail.stdnorminv(t);

  parta = @(nt,a) a^dim .*cos( a*sqrt( nt ));
  partb = @(nt,a) exp(nt*(1-2*a^2)/2);
  fKeister = @(nt,dim,a) parta(nt,a).*partb(nt,a)*(2*pi)^(dim/2);
  y = fKeister(normsqd(yinv(x)),dim,a);
  y(isnan(y)) = 0; % remove any NaN vlaues that could be due to Inf values from stdnorminv
end
