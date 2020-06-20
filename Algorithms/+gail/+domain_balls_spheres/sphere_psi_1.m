function Z = sphere_psi_1(t, d, r, center)
% box-to-sphere transformation
% transformation from a uniform distributions on a (d-1)-dimensional box to a uniform distribution on a d-dimensional sphere.
    t = [zeros(size(t,1),1) t]; %since the first column of the t matrix will be ignored by the TFWW_algorithm function, it is necessary to shift rigth the t matrix
    X = gail.TFWW_algorithm(t,d);% getting points uniformly distributed on a d-dimensional unit sphere 
    Z = bsxfun(@plus, center, r * X); %center + r * X
end