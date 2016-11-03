function [tmu,out_param_AIS]=meanMC_CLT_AIS(Y1,b,d,abstol,alpha,nSig,fudge)
%MEANMC_CLT_AIS combines Monte Carlo method and Adaptive Importance Sampling
%to estimate the mean of a random variable.
%
%   meanMC_CLT_AIS uses adaptive importance sampling to estimate the best
%   value for a transformated variable, i.e., the value that minimizes the
%   variance, and it uses this value to calculate the mean of the function
%   within a specified error tolerance, i.e., | mu - tmu | <= abstol with
%   probability at least 1-alpha, where abstol is the absolute error.
%
%                           Input Arguments
%
%     Y1 --- Anonymous function with two variables, x (independent variable)
%     and b (factor which will be optimized), provided by the user. This
%     function must be the combination between an interest function and
%     the normal density distribution function, i.e, it must be the 
%     importance function.    
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
%     between 0 and 1. Default value is 2e-3.
%
%     alpha --- Uncertainty, which should be a small positive
%     percentage. The default value is 1%.
%
%     nSig --- Number of samples used to compute the sample variance
%
%     fudge --- Standard deviation inflation factor
%
%                           Output Arguments
%
%     tmu --- Estimated mean.
%
%     out_param_AIS.b_value --- value chosen for 'b' within the interval
%     specified.
%
%     out_param_AIS.var --- The sample variance.
%
%     out_param_AIS.sig0 --- standard deviation
%
%     out_param_AIS.ntot --- Total samples used.
%
%     out_param_AIS.time --- The time elapsed in seconds.
%

%     This is a heuristic algorithm based on a Central Limit Theorem
%     approximation combined with adaptive importance sampling.
%
%                             Authors
%
%     BRITO, Rafael de Miranda.PAULO, Ricardo Freitas de.SABARENSE, Mariane de Carvalho.

if nargin < 7
   fudge = 1.2; %variance inflation factor
   if nargin < 6
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 5
         alpha = 0.01; %uncertainty
         if nargin < 4
            abstol = 0.002; %absolute error tolerance
            if nargin < 3
                d = 1; %number of dimensions
                if nargin < 2
                    b=[-2 2]; % 'b'interval
                    if nargin < 1 
                        Y1=input('Please inform "g(x)"');% 'g(x)'
                    end                  
                end
            end
         end
      end
   end
end
  
%                               RESTRICTIONS

% Checking function input

if isa(Y1,'function_handle') == 0 || nargin(Y1) ~= 2
    warning('meanMC_AIS_g:Y1notafunction',...
    ['"Y" must be a function handle with two variables - "n" and "b".\n'...
        'A default function "Y(n,b) = (n+b).^2" will be used']);
    Y1 =@(n,b)(n+b).^2;
end

% Checking 'b' input

if isa(b,'double') == 0 || numel(b) ~= 2 || issorted(b) == 0
    warning('meanMC_AIS_g:invalidInterval',...
    ['"b" must be an array with 2 elements in ascending order.\n'...
    'A default interval of [-2 2] will be used.']);
    b=[-2 2];
    
end

% Checking dimension input

if d <= 0 || mod(d,1) ~= 0 
    warning('meanMC_AIS_g:invalidDimension',...
    ['The number of dimensions must be an integer positive number.\n'...
    'A default value d=1 will be used.']);
    d=1;
end

% Checking other inputs

if abstol <= 0 || abstol >1
    warning('meanMC_AIS_g:invalidTolerance',...
        ['The absolut tolerance must be between zero and one.\n'...
    'A default value abstol = 0.002 will be used.']);
    abstol = 0.002;
end

if alpha <= 0
    warning('meanMC_AIS_g:invalidUncertainty',...
        ['The uncertainty must be higher than zero.\n'...
    'A default value alpha = 0.01 will be used.']);
    alpha = 0.01;
end

if nSig <= 0
    warning('meanMC_AIS_g:invalidNumberSamples',...
        ['The number of samples must be higher than zero.\n'...
    'A default value nSig = 1e4 will be used.']);
    nSig = 1e4;
end

if fudge <= 0
    warning('meanMC_AIS_g:invalidInflationFactor',...
        ['The inflation factor must be higher than zero\.'...
    'A default value fudge = 1.2 will be used.']);
    fudge = 1.2;
end
%__________________________________________________________________________

out_param_AIS.alpha = alpha; % Save the input parameters to a structure.
out_param_AIS.abstol = abstol;
out_param_AIS.fudge = fudge;
out_param_AIS.nSig = nSig;

b_vec=[b(1),((b(1)+b(2))/2),b(2)]; % Generates a vector with 3 values equally spaced
%within the interval defined.

%Y = @(n,b)Y1(randn(n,d),b); % Integrand evaluated at the sample points. 

tstart=tic; % Starts the clock.

var_b=b_vec; % Avoids the change of vector size inside the loop.

% Checking the variance for each element in b_vec:
for i=1:numel(b_vec)
    b_value=b_vec(i);    
    var_b(i)=var(Y1(randn(nSig,d),b_value));
end
[S_var,S_pos]=min(var_b); % Computes the position and the value of the smaller
% variance calculated from 'b_vec'

% Parabolic interpolation between b_vec and calculated variance:
A=[b_vec'.^2 b_vec' ones(3,1)];
p=A\var_b';

% Minimum search using the approximated parabola:
[x]=(-p(2)/(2*p(1)));
% Variance calculation using the value estimated as the minimum 'x'
var_bx=var(Y1(nSig,x)); 


% Checking the best value for b:
if var_bx < S_var && var_bx > 0
    out_param_AIS.b_value = x;
    out_param_AIS.var = var_bx;
else
    out_param_AIS.b_value = b_vec(S_pos);
    out_param_AIS.var = S_var;
end

out_param_AIS.sig0 = sqrt(out_param_AIS.var); %standard deviation
sig0up = out_param_AIS.fudge.*out_param_AIS.sig0; %upper bound on the standard deviation
nmu = max(1,ceil((-norminv(alpha)*sig0up/abstol).^2)); 
   %number of samples needed for mean
assert(nmu<1e8,['nmu = ' int2str(nmu) ', which is too big']) 
   %don't exceed sample budget
tmu=mean(Y1(randn(nmu,1),out_param_AIS.b_value)); %estimated mean
out_param_AIS.ntot=4.*nSig+nmu; %total samples required
out_param_AIS.time=toc(tstart); %elapsed time    

end

