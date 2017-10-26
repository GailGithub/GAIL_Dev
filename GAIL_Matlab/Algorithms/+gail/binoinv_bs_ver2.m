function x = binoinv_bs_ver2(y,n,p)
%BINOINV_BS Inverse of the binomial CDF using binary search
%   X = BINOINV_BS(Y,N,P) returns the inverse of the binomial cdf with 
%   parameters N and P. Since the binomial distribution is
%   discrete, BINOINV returns the least integer X such that 
%   the binomial cdf evaluated at X, equals or exceeds Y.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Note that X takes the values 0,1,2,...,N.
%
%   See also BINOINV, BINOCDF, BINOFIT, BINOPDF, BINORND, BINOSTAT, ICDF.

% Nade Sritanyaratana
% Created: August 28, 2014
% v1.0 - Obtained documentation from BINOINV.

x = zeros(size(y)); %initialize vector of answers
for i=1:numel(y)
    % Initialize values for the current iteration
    yi = y(i); %value of probability to match
    xl = 0; % answer is between xl
    xr = n; % and xr;
    yl = binocdf(xl, n, p); %cdf at the left endpoint is correct
    ni = xr-xl; %width of the interval
    if yl-yi>=0 %is the left endpoint is correct
       x(i) = xl; %assign answer
    else
       % Binary search      
       while ni>1 % xl and xr not adjacent
           xmid = xl+floor(ni/2); %calculate midpoint
           ymid = binocdf(xmid, n, p); %new cdf falue
           if ymid-yi >= 0 %midpoint is the new right
              xr=xmid;
           else %midpoint is the new left
              xl=xmid;
           end
           ni = xr-xl; %new interval width
       end
       x(i) = xr; %assign answer
    end  
end
end