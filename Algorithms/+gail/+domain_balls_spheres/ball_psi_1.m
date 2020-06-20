function Z = ball_psi_1(t, d, r, center)
% box-to-ball transformation function
% transformation from a uniform distributions on a d-dimensional box to a uniform distribution on a d-dimensional ball.
    X = gail.TFWW_algorithm(t,d);% getting points uniformly distributed on a d-dimensional sphere 
    u = t(:,1);% a uniform distribution on the interval [0, 1]
    % now, a sphere-to-ball transformation must be applied
    Z = bsxfun(@plus, center, r * bsxfun(@times, u.^(1.0/d), X)); %center + r * ((u^(1/d))*t)
end