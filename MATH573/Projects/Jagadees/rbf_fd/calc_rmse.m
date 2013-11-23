% CALC_RMSE calculates Root mean square error for the approximated function
%         with respect to the original function
%
%    [RMSE, RMSED, STABLE_FLAG, U_HAT] = CALC_RMSE(DM, POINTS, INT_PTS, BOUND_PTS, 
%                                     G, U, LU, EPSILON, RBFQR_FLAG)
%    Input:
%    dm : is distance matrix which lists euclidean distance between point(i)
%         to point(j). 
%    points: is the index to the set of points to be used for solving
%               the equation
%    int_pts: index to the set of interior points
%    bound_pts: index to the set of boundary points
%    g: function defining the boundary behaviour
%    u: actual function
%    Lu: Laplacian of the function 'u' to be approximated
%    epsilon: shape parameter for Gaussian RBF
%    RBFQR_flag: Flag if 'true' allows using stable computation else only
%                direct computation used
%    Output:
%    RMSE: RMS Error for function approximation
%    RMSED: RMS Error for differentiation
%    Stable_Flag = Indicating if the stable computation used for Weight computing
%    u_hat: Approximated function values

function [rmse, rmsed, stable_flag, u_hat] = calc_rmse(dm, points, int_pts, bound_pts, g, u, Lu, epsilon, RBFQR_flag)

    % Compute A Matrix to find the approximation to function 'u'
    [A, b, stable_flag] = solve_poisson(dm, points, int_pts, bound_pts, g, Lu, epsilon, RBFQR_flag);

    %fprintf('cond(A) %f\n, cond(A));
    u_hat = (A\b);
    Lu_exact = zeros(length(points), 1);
    u_exact  =  u(points(1,:), points(2,:))';
    Lu_exact(:) = Lu(points(1,:), points(2,:))';
    Lu_hat = A*u_exact;

    rmse  = norm(u_exact(int_pts)-(u_hat(int_pts)))/sqrt(length(int_pts));
    rmsed = norm(Lu_exact(int_pts)-(Lu_hat(int_pts)))/sqrt(length(int_pts));
end

