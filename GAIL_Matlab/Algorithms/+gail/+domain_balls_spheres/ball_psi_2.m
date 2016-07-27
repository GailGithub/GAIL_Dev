function Z = ball_psi_2(t, d, r, center)
% normal-to-ball transformation
% transformation from a normal distributions on a d-dimensional space to a uniform distribution on a d-dimensional ball.
    radius2 = @(t) sum(t.^2, 2); %computes the radius squared
    factor = @(t) ((chi2cdf(radius2(t),d).^(1.0/d))./sqrt(radius2(t)));% computes the factor by which the input must be multiplied

    Z = bsxfun(@plus, center, r * bsxfun(@times,t,factor(t)));% center + r * (factor*t)
end