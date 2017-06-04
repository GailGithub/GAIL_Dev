% productfun produces the product function used to section 3.5.1 in Lan
% Jiang's thesis.
function output = productfun(center,hyperbox,inparam)

[output.iidQ, output.iidparam] = ...
    cubMC_g(@(x)productFun(x,center),hyperbox,inparam);
[output.SobolQ, output.Sobolparam] = ...
    cubSobol_g(@(x)productFun(x,center) , hyperbox, inparam);
[output.LatticeQ, output.Latticeparam] =...
    cubLattice_g(@(x)productFun(x,center), hyperbox, inparam);
 
end
function f_eval = productFun(x,center)
f_eval = ones(size(x,1),1);
for i = 1:size(center,2)
f_eval = f_eval.*(x(:,i).^2+center(i));
end
end
 
% function f_eval = transMVNP(d,e,box,w,C)
%     f_eval = (e-d)*ones(size(w,1),1);
%     aux = ones(size(w,1),1);
%     y = [];
%     for i = 2:size(box,2);
%         y = [y gail.stdnorminv(d+w(:,i-1).*(e-d))];
%         aux = sum(bsxfun(@times,C(i,1:i-1),y),2);
%         a = (box(1,i)-aux)/C(i,i);
%         b = (box(2,i)-aux)/C(i,i);
%         d = gail.stdnormcdf(a);
%         e = gail.stdnormcdf(b);
%         f_eval = f_eval .* (e-d);
%     end
% end