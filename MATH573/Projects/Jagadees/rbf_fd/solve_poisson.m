% SOLVE_POISSON used for approximation of the Dirichlet problem for the Poisson equation 
%               It Gaussian radial basis function finite difference (RBF-FD)
%               method with irregular centres.
%
%   [A, b] = SOLVE_POISSON(dm, points, int_pts, bound_pts, g, Lu, epsilon, 
%            RBFQR_flag) return the Stencil support where
%    dm : is distance matrix which lists euclidean distance between point(i)
%         to point(j). 
%    points: is the index to the set of points to be used for solving
%               the equation
%    g: function defining the boundary behaviour
%    Lu: Laplacian of the function 'u' to be approximated
%    int_pts: index to the set of interior points
%    bound_pts: index to the set of boundary points
%    epsilon: shape parameter for Gaussian RBF
%    RBFQR_flag: Flag if 'true' allows using stable computation else only
%                direct computation used
%
%   See also FIND_STENCIL_WEIGHTS
%
% Reference: O. Davydov, D.T. Oanh, On the optimal shape parameter for 
% Gaussian radial basis function finite difference approximation of the 
% Poisson equation, Comput. Math. Appl. 62 (2011) 2143–2161.

function [A, b, stable_flag] = solve_poisson(dm, points, int_pts, bound_pts, g, Lu, epsilon, RBFQR_flag)

    A = zeros(length(points'), length(points'));
    b = zeros(length(points'), 1);
    v = zeros(length(points'), 1);
    stable_flag = false;
    
    for i=bound_pts
        A(i,i) = 1.0;
    end

    b(int_pts)   = Lu(points(1,int_pts), points(2,int_pts));
    b(bound_pts) = g(points(1,bound_pts), points(2,bound_pts));
    
    % Computes Stencil points and then Stencil weights to populate the A
    % Matrix to solve the equation
    for i=int_pts
        center = i; %index of the point to be used as center
        
        % Find the Stencil points
        stencil_support = stencil_support_selection(dm, points, center);

        % compute Stencil weights
        [weights, v(i), stable_flag_i] = find_stencil_weights(points(:,stencil_support)', epsilon, RBFQR_flag);
        stable_flag = stable_flag || stable_flag_i;
        A(i,stencil_support) = real(weights);
    end
end
