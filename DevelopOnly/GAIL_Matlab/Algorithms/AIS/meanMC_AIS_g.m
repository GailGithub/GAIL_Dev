function [tmu,out_param]=meanMC_AIS_g(Y,b_vec,d,abstol,alpha,nSig,fudge)

%   MEANMC_AIS_g uses adaptive importance sampling to estimate the value of
%   a function integral within a specified error tolerance, i.e., 
%   | mu - tmu |<= abstol with probability at least 1-alpha, where abstol
%   is the absolute error tolerance using the adaptive importance sampling.  
%
%
%                           Input Arguments
%
%     Y --- Anonymous function provided by the user.
%     
%     B --- Vector with three values for standard deviation.
%
%     d --- number of dimensions.
%
%     abstol --- Absolute error tolerance, which should be
%     positive, default value is 2e-3.
%
%     alpha --- Uncertainty, which should be a small positive
%     percentage. The default value is 1%.
%
%     nSig --- Number of samples used to compute the sample variance
%
%     fudge --- Standard deviation inflation factor
%
%
%                           Output Arguments
%
%
%     tmu --- the estimated value of the integral.
%
%     out_param.ntot --- total samples used.
%
%     out_param.var --- the sample variance.
%
%     out_param.time --- the time elapsed in seconds.
%



if nargin < 6
   fudge = 1.2; %variance inflation factor
   if nargin < 5
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 4
         alpha = 0.01; %uncertainty
         if nargin < 3
            abstol = 0.002; %absolute error tolerance
            if nargin < 3
                d = 3; %number of dimensions
                if nargin < 2
                    b_vec=[0.5,1.5,2.5]; 
                    if nargin < 1
                        Y1 = @(x,b) ((sqrt(2.*pi).*b).^d).*cos(b.*sqrt(sum(x.*x,2)))...
                            .*exp((1/2-b.^2).*sum(x.*x,2)); %integrand
                        Y = @(n,b)Y1(randn(n,d),b); % integrand evaluated at
                        % the sample points. 
                    end                  
                end
            end
         end
      end
   end
end
  
% RESTRICTIONS

if isa(Y,'function_handle') == 0
    error('"Y" must be a function handle');
end

if isa(b_vec,'double') == 0
    error('"b_vector" must be an array');
end

if d <= 0 
    error('The number of dimensions must be an integer positive number');
end


out_param.alpha = alpha; %save the input parameters to a structure
out_param.fudge = fudge;
out_param.nSig = nSig;

tstart=tic; %starts the clock.

var_b=b_vec;
for i=1:numel(b_vec)
    b_value=b_vec(i);    
    var_b(i)=var(Y(nSig,b_value));
end

%parabolic interpolation:
A=[b_vec'.^2 b_vec' ones(3,1)];
p=A\var_b';
fmin=@(x)p(3)+p(2)*x+p(1)*(x.^2);
[x]=fminbnd(fmin,b_vec(1),b_vec(3));
var_bx=var(Y(nSig,x));

[S_var,S_pos]=min(var_b);

%Checking the best value for b:
if var_bx < S_var
    out_param.b_value = x;
    out_param.var = var_bx;
else
    out_param.b_value = b_vec(S_pos);
    out_param.var = S_var;
end

% MeanMC_g calculation



tmu = meanMC_g(@(n)Y(n,out_param.b_value),abstol,0);
    

sig0 = sqrt(out_param.var); %standard deviation
out_param.time=toc(tstart); %elapsed time    



end

