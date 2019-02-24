
function output = MultivarNorProb(box,mu,cov,inparam)
% MultivarNorProb computes the cumulative distribution function of the
% multivariate normal with mean mu, covariance matrix cov and within the
% region defined by box.
box = bsxfun(@minus, box,mu');
C = chol(cov)'; s = size(C,1);
a = box(1,1)/C(1,1); 
b = box(2,1)/C(1,1);
d = gail.stdnormcdf(a); 
e = gail.stdnormcdf(b);
transhyperbox = [zeros(1,s-1);ones(1,s-1)];
[output.iidQ, output.iidparam] = ...
    cubMC_g(@(x) transMVNP(d,e,box,x,C),transhyperbox,inparam);
[output.SobolQ, output.Sobolparam] = ...
    cubSobol_g(@(x) transMVNP(d,e,box,x,C) , transhyperbox, inparam);
[output.LatticeQ, output.Latticeparam] =...
    cubLattice_g(@(x) transMVNP(d,e,box,x,C), transhyperbox, inparam);
 
end
 
function f_eval = transMVNP(d,e,box,w,C)
    f_eval = (e-d)*ones(size(w,1),1);
    aux = ones(size(w,1),1);
    y = [];
    for i = 2:size(box,2);
        y = [y gail.stdnorminv(d+w(:,i-1).*(e-d))];
        aux = sum(bsxfun(@times,C(i,1:i-1),y),2);
        a = (box(1,i)-aux)/C(i,i);
        b = (box(2,i)-aux)/C(i,i);
        d = gail.stdnormcdf(a);
        e = gail.stdnormcdf(b);
        f_eval = f_eval .* (e-d);
    end
end