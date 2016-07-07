function X = TFWW_algorithm(t, d)
% algorithm that computes the transformation from a uniform distributions on a d-dimensional box
% to a uniform distribution on a d-dimensional sphere, using the d-1 last coordinates of t
    if mod(d,2) == 0 % d even
        m = d/2;
        g = zeros(size(t,1),m+1);% g is 1-based, diferently from the original algorithm
        g(:,m+1) = ones(size(t,1),1);
        
        %computing g(:,i) = g(:,i+1).*(t(:,i).^(1.0/(i-1))) without for
        %loops for efficiency
        exponents = ones(1, m-1)./(1:m-1);
        g(:,2:m) = cumprod(bsxfun(@power, t(:,2:m), exponents), 2, 'reverse');
        
        %computing dd(:,i) = sqrt(g(:,i+1) - g(:,i)) without for loops
        dd = sqrt(diff(g,1,2));
        
        %computing the following without for loops:
        %   X(:,2*i-1) = dd(:,i).*cos(2*pi*t(:,m+i));
        %   X(:,2*i) = dd(:,i).*sin(2*pi*t(:,m+i));
        X = reshape([dd.*cos(2*pi*t(:,m+1:2*m)); dd.*sin(2*pi*t(:,m+1:2*m))], size(t,1), 2*m);

    else %d odd
        m = (d-1)/2;
        g = zeros(size(t,1),m+1);% g is 1-based, diferently from the original algorithm
        g(:,m+1) = ones(size(t,1),1);
        
        %computing g(:,i) = g(:,i+1).*(t(:,i).^(2.0/(2*i-1))) without for
        %loops for efficiency
        exponents = (2*ones(1, m-1))./(3:2:2*m-1);
        g(:,2:m) = cumprod(bsxfun(@power, t(:,2:m), exponents), 2, 'reverse');
        
        %computing dd(:,i) = sqrt(g(:,i+1) - g(:,i)) without for loops
        dd = sqrt(diff(g,1,2)); % dd is equivalent to the d vector in the original algorithm
              
        X = zeros(size(t)); %initializing the answer matrix for efficiency 
        
        X(:,1) = dd(:,1).*(1-2*t(:,m+1));
        X(:,2) = 2*dd(:,1).*sqrt(t(:,m+1).*(1-t(:,m+1))).*cos(2*pi*t(:,m+2));%multiplied by 2 to match formula (cf.(1.5.28))
        X(:,3) = 2*dd(:,1).*sqrt(t(:,m+1).*(1-t(:,m+1))).*sin(2*pi*t(:,m+2));%multiplied by 2 to match formula (cf.(1.5.28))
        
        %computing the following without for loops:
        %   X(:,2*i) = dd(:,i).*cos(2*pi*t(:,m+i+1));
        %   X(:,2*i+1) = dd(:,i).*sin(2*pi*t(:,m+i+1));
        if d > 3
            X(:,4:end) = reshape([dd(:,2:end).*cos(2*pi*t(:,m+3:2*m+1)); dd(:,2:end).*sin(2*pi*t(:,m+3:2*m+1))], size(t,1), 2*(m-1));
        end
    end
end
