function [tmu,out_param]=meanMC_AIS_g(Y1,b,d,abstol,alpha,nSig,fudge)

%   MEANMC_AIS_g uses adaptive importance sampling to estimate the value of
%   a function integral within a specified error tolerance, i.e., 
%   | mu - tmu |<= abstol with probability at least 1-alpha, where abstol
%   is the absolute error tolerance using the adaptive importance sampling.  
%
%
%                           Input Arguments
%
%     Y --- Anonymous function with two variables, x (independent variable)
%     and b (factor which will be optimized), provided by the user. This
%     function must be the combination between an interest function and
%     the normal density distribution function - the importance function.
%     
%     b --- Vector with two values that indicate an interval to be used for
%     variable 'b'. This interval will be used to create a vector, b_vec,
%     containing three equally spaced points which will be used to generate 
%     a parabolic interpolation (between b_vec and the corresponding 
%     variance) to determine a local minimum within the interval specified.
%
%     d --- Number of dimensions.
%
%     abstol --- Absolute error tolerance, which should be
%     positive, default value is 2e-3.
%
%     alpha --- Uncertainty, which should be a small positive
%     percentage. The default value is 1%.
%
%     nSig --- Number of samples used to compute the sample variance.
%
%     fudge --- Standard deviation inflation factor.
%
%
%                           Output Arguments
%
%
%     tmu --- Estimated value of the integral.
%
%     out_param.ntot --- Total samples used.
%
%     out_param.var --- Variance.
%
%     out_param.time --- Time elapsed (in seconds).
%
%     Example
%
%    Calculate the mean of 
%
%                             Authors
%
%     BRITO, Rafael de Miranda.
%     Ricardo.
%     SABARENSE, Mariane de Carvalho.


if nargin < 6
   fudge = 1.2; %variance inflation factor
   if nargin < 5
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 4
         alpha = 0.01; %uncertainty
         if nargin < 3
            abstol = 0.002; %absolute error tolerance
            if nargin < 3
                d = 1; %number of dimensions
                if nargin < 2
                    b=[0.5 1.5]; 
                    if nargin < 1
                        Y1 = @(x,b) ((sqrt(2.*pi).*b).^d).*cos(b.*sqrt(sum(x.*x,2)))...
                            .*exp((1/2-b.^2).*sum(x.*x,2)); % g(x)
                    end                  
                end
            end
         end
      end
   end
end
  
%                               RESTRICTIONS

% Checks function input

if isa(Y1,'function_handle') == 0 || nargin(Y1) ~= 2
    error('"Y" must be a function handle with two variables - "x" and "b".');
end

% Checks 'b' input

if isa(b,'double') == 0 || numel(b) ~= 2 || issorted(b) == 0
    error('"b" must be an array with 2 elements in ascending order');
    
end

% Checks dimension input

if d <= 0 || mod(d,1) ~= 0 
    error('The number of dimensions must be an integer positive number');
end
%__________________________________________________________________________


out_param.alpha = alpha; % Save the input parameters to a structure.
out_param.fudge = fudge;
out_param.nSig = nSig;

b_vec=linspace(b(1),b(2),3); % Generates a vector with 3 values equally spaced
%   within the interval defined.

Y = @(n,b)Y1(randn(n,d),b); % g(x) evaluated at the sample points. 

tstart=tic; % Starts the clock.

var_b=b_vec; % Avoids the change of vector size inside the loop.

% Checking the variance for each element in b_vec:
for i=1:numel(b_vec)
    b_value=b_vec(i);    
    var_b(i)=var(Y(nSig,b_value));
end
[S_var,S_pos]=min(var_b); % Computes the position and the value of the smaller
% variance calculated from 'b_vec'

% Parabolic interpolation between b_vec and calculated variance:
A=[b_vec'.^2 b_vec' ones(3,1)];
p=A\var_b';
fmin=@(x)p(3)+p(2)*x+p(1)*(x.^2);

% Minimum search using the approximated parabola:
[x]=fminbnd(fmin,b_vec(1),b_vec(3));

% Variance calculation using the value estimated as the minimum 'x'
var_bx=var(Y(nSig,x)); 


% Checking the best value for b:
if var_bx < S_var
    out_param.b_value = x;
    out_param.var = var_bx;
else
    out_param.b_value = b_vec(S_pos);
    out_param.var = S_var;
end

% MeanMC_g calculation

[tmu,out_param]=meanMC_g(@(n)Y(n,out_param.b_value),abstol,0);

out_param.nTotal= 4.*nSig+(out_param.ntot);
sig0 = sqrt(out_param.var); %standard deviation
out_param.time=toc(tstart); %elapsed time    



end

