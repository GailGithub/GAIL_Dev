
function fval = GenzFunc(w,params)
dim = numel(params.a);
nn = size(w,1);
am = params.a - params.mu;
bm = params.b - params.mu;
a1 = am(1)/params.CovProp.C(1,1);
b1 = bm(1)/params.CovProp.C(1,1);
d = gail.stdnormcdf(a1);
e = gail.stdnormcdf(b1);
fval = (e-d)*ones(nn,1);
y = zeros(nn,dim-1);
for i = 2:dim
  y(:,i-1) = gail.stdnorminv(d+w(:,i-1).*(e-d));
  aux = sum(bsxfun(@times,params.CovProp.C(i,1:i-1),y(:,1:i-1)),2);
  a1 = (am(i)-aux)/params.CovProp.C(i,i);
  b1 = (bm(i)-aux)/params.CovProp.C(i,i);
  d = gail.stdnormcdf(a1);
  e = gail.stdnormcdf(b1);
  fval = fval .* (e-d);
end

fval(isnan(fval)) = 0; % reset NanN vlaues to zero
end
