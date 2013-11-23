% STENCIL_SUPPORT_SELECTION used to find the 2D Stencil support
%   which will be used for Finite difference Approxmn
%
%   Set_Xi_zeta = STENCIL_SUPPORT_SELECTION(dm, point_set, center,
%   bound_pts) return the Stencil support where
%    dm : is distance matrix which lists euclidean distance between point(i)
%         to point(j). 
%    point_set: is the index to the set of points to be used for stenil
%               selection
%    center: index to the center point
%    bound_pts: index to the set of boundary points
%    Set_Xi_zeta: index to the set of stencil points
%
%   See also FIND_STENCIL_WEIGHTS
%
% Reference: O. Davydov, D.T. Oanh, On the optimal shape parameter for 
% Gaussian radial basis function finite difference approximation of the 
% Poisson equation, Comput. Math. Appl. 62 (2011) 2143–2161.

%Algorithm 1. Meshless stencil support selection
function Set_Xi_zeta = stencil_support_selection(dm, point_set, center)

	k = 6;
	m = 30;
	u = 3.0;

	% I. Find m nearest points
	m_near_pts_I = get_near_points(dm, center, m);
    m_near_points_polar = conv2polar(point_set, m_near_pts_I, center);
    
 	Set_Xi_zeta = 1:k; %Initialise the set to k nearest pts
    alpha = get_alpha(Set_Xi_zeta, m_near_points_polar);
    if (max(alpha) <= u*min(alpha))
        Set_Xi_zeta = [center, m_near_pts_I(Set_Xi_zeta)];
        return;
    end
   
    cost = dot(alpha, alpha);
    %fprintf('c: %d, Initial cost %f\n', center, cost)
    
    % II. For i=n+1,...,m:
    for i=k+1:m

        % 1. Extended Set
        Set_Xi_zeta_extend = [Set_Xi_zeta, i];
        [alpha, kI] = get_alpha(Set_Xi_zeta_extend, m_near_points_polar);
        
        % 2.
        ii = find(kI==i);
        if ii==1
            ai1 = alpha(end);
        else
            ai1 = alpha(ii-1); %angle between ith and prior neighbor (counter clockwise)
        end
        ai2 = alpha(ii);   %angle between ith and next neighbor (counter clockwise)
        
        min_alpha = min(alpha);
        if ai1 > min_alpha && ai2 > min_alpha

            % find the pth elem to be removed
            p = find_pth_elem2remv(alpha);

            elem_to_remv = kI(p); %get the point to remove
            
            % ii check the cost of new stencil
            Set_Xi_zeta_new = setdiff(Set_Xi_zeta_extend, elem_to_remv);
            alpha1 = get_alpha(Set_Xi_zeta_new, m_near_points_polar);
            new_cost = dot(alpha1, alpha1);
            if new_cost < cost
                %fprintf('new_cost %f\n', new_cost);
                Set_Xi_zeta = Set_Xi_zeta_new;
                cost = new_cost;
                if (max(alpha1) <= u*min(alpha1))
                    Set_Xi_zeta = [center, m_near_pts_I(Set_Xi_zeta)];
                    return;
                end
            else
                %fprintf('New cost is not lesser: check\n');
            end
            
        end

    end
    
    % Algorithm has not terminated yet, in the previous loop
    % III
    [alpha, kI] = get_alpha(Set_Xi_zeta, m_near_points_polar);
    % Here (max(alpha) > u*min(alpha))
    
    % find the pth elem to be removed
    p = find_pth_elem2remv(alpha);
    
    elem_to_remv = kI(p); %get the point to remove
    Set_Xi_zeta_new = setdiff(Set_Xi_zeta, elem_to_remv);
    Set_Xi_zeta = [center, m_near_pts_I(Set_Xi_zeta_new)];
    fprintf('Algo fails to find the optimal stencil, Default returned !!\n');
    return;
    
end

% Find the the point to remove which improves the angle spreading in the
% stencil
function p = find_pth_elem2remv(alpha)
    j = find(alpha == min(alpha));
    if length(j) > 1
        j = j(1);
    end
    % find the pth elem to be removed
    if j==1
        a_jn1 = alpha(end);
        a_jp1 = alpha(j+1);
    elseif j==length(alpha)
        a_jn1 = alpha(j-1);
        a_jp1 = alpha(1);
    else
        a_jn1 = alpha(j-1);
        a_jp1 = alpha(j+1);
    end
    if a_jn1 < a_jp1
        p = j;
    else
        p = j+1;
        if j==length(alpha)
            p = 1;
        end
    end
end

function m_near_points_polar = conv2polar(point_set, m_near_pts_I, center)
    % Transform all the points wrt center
    % convert to polar coordinates
    m_near_points_polar = zeros(size(point_set(:,m_near_pts_I)));
    for i=1:length(m_near_pts_I)
        xy = point_set(:,m_near_pts_I(i)) - point_set(:,center);
        %[THETA,RHO] = cart2pol(X,Y)
        [theta, rho] = cart2pol(xy(1), xy(2));
        theta = theta*180/pi;
        if theta<0
            theta = theta+360;
        end
        m_near_points_polar(:,i) = [theta;rho];
    end
end

% finds angle between each two neighboring points, angle bearing at the
% center
function [alpha, kI] = get_alpha(Set_Xi_zeta, m_near_points_polar)

    thetas = m_near_points_polar(1,Set_Xi_zeta);
    [t,I] = sort(thetas);
    alpha = zeros(length(thetas), 1);
    
    for i=1:length(thetas)
        ip1 = 1 + mod(i, length(thetas));
        alpha(i) = t(ip1)-t(i);
        if alpha(i) < 0
            alpha(i) = 360+alpha(i);
        end
    end
    if abs(sum(alpha) - 360.0) > 1e-3
        error('Wrong Alpha computation, should sum to 360')
    end
    kI = Set_Xi_zeta(I);

end

% finds the distance from the center to the each stencil points
function m_near_pts = get_near_points(dm, center, m)

    dmc = dm(center,:);
    [B,I] = sort(dmc);
    assert(B(1)==0);
    m_near_pts = I(2:m+1);
    
end
