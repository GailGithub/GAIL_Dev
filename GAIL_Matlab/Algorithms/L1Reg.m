% L1LinearRegression: Calculates L-1 multiple linear regression by IRLS
% by Will Dwinnell
%
% x = L1LinearRegression(A,b)
%
% x  = discovered linear coefficients
% A  = independent variables
% b  = dependent variable
%
%
% Last modified: Mar-27-2009

function x = L1LinearRegression(A,b)

% Determine size of predictor data
[n m] = size(A);

% Initialize with least-squares fit
x       = A\ b;  % Least squares regression
xOld    = x;
xOld(1) = xOld(1) + 1e-4;  % Force divergence

% Repeat until convergence
while (max(abs(x - xOld)) > 1e-8)
    % Move old coefficients
    xOld = x;
    
    % Calculate new observation weights (based on residuals from old coefficients)
    W = sqrt(1 ./ max(abs(( A * xOld) - b),1e-8));  % Floor to avoid division by zero
    
    % Calculate new coefficients
    x = (repmat(W,[1 m]) .* A) \ (W .* b);
end


% EOF


