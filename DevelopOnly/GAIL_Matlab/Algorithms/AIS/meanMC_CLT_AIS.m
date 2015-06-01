function [tmu,out_param]=meanMC_CLT_AIS(fx,abstol,alpha,nSig,fudge)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   tmu = MEANMC_CLT(Yrand,abstol,alpha,nSig,fudge) estimates the mean, mu, of a random variable Y to
%   within a specified error tolerance, i.e., | mu - tmu | <= abstol with
%   probability at least 1-alpha, where abstol is the absolute error
%   tolerance.  The default values are abstol=1e-2 and alpha=1%. Input
%   Yrand is a function handle that accepts a positive integer input n and
%   returns an n x 1 vector of IID instances of the random variable Y.
%
%   Input Arguments
%
%     Yrand --- the function for generating n IID instances of a random
%     variable Y whose mean we want to estimate. Y is often defined as a
%     function of some random variable X with a simple distribution. The
%     input of Yrand should be the number of random variables n, the output
%     of Yrand should be n function values. For example, if Y = X.^2 where X
%     is a standard uniform random variable, then one may define Yrand =
%     @(n) rand(n,1).^2.
%     
%     fx --- defined function
%
%     abstol --- the absolute error tolerance, which should be
%     positive, default value is 1e-2.
%
%     alpha --- the uncertainty, which should be a small positive
%     percentage. The default value is 1%.
%
%     nSig --- the number of samples used to compute the sample variance
%
%     fudge --- the standard deviation invlation factor
%
%   Output Arguments
%
%     tmu --- the estimated mean of Y.
%
%     out_param.ntot --- total sample used.
%
%     out_param.var --- the sample variance.
%
%     out_param.time --- the time elapsed in seconds.
%

%This is a heuristic algorithm based on a Central Limit Theorem
%approximation
if nargin < 5
   fudge = 1.2; %variance inflation factor
   if nargin < 4;
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 3
         alpha = 0.01; %uncertainty
         if nargin < 2
            abstol = 0.01; %absolute error tolerance
            if nargin < 1
                fx=@(x)max(x,0);
            end
         end
      end
   end
end
nMax=1e8; %maximum number of samples allowed.
out_param.alpha = alpha; %save the input parameters to a structure
out_param.fudge = fudge;
out_param.nSig = nSig;

b=[-1,0,1];
var_b=b;
%dens_func = @(b)(1/(sqrt(2.*pi)).*exp(-(b.^2)/2));
%normal distribution density function
tstart=tic; %start the clock
fx_b=@(t,b_value)fx(t+b_value).*exp(-t.*b_value-b_value.^2/2);
for i=1:numel(b)
    b_value=b(i);
    var_b(i)=var(fx_b(randn(nSig,1),b_value));
end

%d_value=dens_func(b);
    %Yval = x_values(nSig);% get samples to estimate variance 
    %fx=max(Yval,0);
    %integrand=fx*d_value;
    %u=mean(integrand);
    %variance=var(integrand);
[S_var,S_pos]=min(var_b);
out_param.var = S_var; %calculate the sample variance--stage 1
out_param.b_value=b(S_pos);
sig0 = sqrt(out_param.var); %standard deviation
sig0up = out_param.fudge.*sig0; %upper bound on the standard deviation
nmu = max(1,ceil((-norminv(alpha)*sig0up/abstol).^2)); 
   %number of samples needed for mean
assert(nmu<nMax,['nmu = ' int2str(nmu) ', which is too big']) 
   %don't exceed sample budget
tmu=mean(fx_b(randn(nmu,1),out_param.b_value)); %estimated mean
out_param.ntot=nSig+nmu; %total samples required
out_param.time=toc(tstart); %elapsed time    


%The real value can be calculated through the integral of the following
%expression between the interval [-3 3]:

%f(x)*p(x)=x.*(1./(2.*pi)).*exp(-(x.^2)/2)

% integral = (-exp(-(x^2)/2))/(2*pi)

%The result for the interval specified is equal to zero.


%out_param.var = var(integrand); %calculate the sample variance--stage 1
%sig0 = sqrt(out_param.var); %standard deviation
%sig0up = out_param.fudge.*sig0; %upper bound on the standard deviation
%nmu = max(1,ceil((-norminv(alpha)*sig0up/abstol).^2)); 
   %number of samples needed for mean
%assert(nmu<nMax,['nmu = ' int2str(nmu) ', which is too big']) 
   %don't exceed sample budget
%tmu=mean(Yrand(nmu)); %estimated mean
%out_param.ntot=nSig+nmu; %total samples required
toc(tstart); %elapsed time
end

