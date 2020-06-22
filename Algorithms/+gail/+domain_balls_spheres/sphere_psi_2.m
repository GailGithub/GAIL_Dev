function Z = sphere_psi_2(t, d, r, center)
% normal-to-sphere transformation
% transformation from a normal distributions on a d-dimensional space to a uniform distribution on a d-dimensional sphere.
    radius = sqrt(sum(t.^2, 2)); %computes the radius
    Z = bsxfun(@plus, center, r * bsxfun(@times,t,radius.^(-1)));% center + r * (factor*t)
end