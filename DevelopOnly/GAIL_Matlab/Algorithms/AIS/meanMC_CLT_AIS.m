function [tmu,out_param]=meanMC_CLT_AIS(fx,abstol,alpha,nSig,fudge)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   tmu = meanMC_CLT(fx,abstol,alpha,nSig,fudge) estimates the mean, of a random fuction to
%   within a specified error tolerance, i.e., | mu - tmu | <= abstol with
%   probability at least 1-alpha, where abstol is the absolute error
%   tolerance using the adaptive importance sampling.  
%   The default values are abstol=1e-2 and alpha=1%. Input 'fx' 
%   is a function handle that accepts a positive integer input n and
%   returns an n x 1 vector of pseudorandom values drawn
%   from the standard normal distribution of the random function.
%
%                           Input Arguments
%
%     fx --- the function for generating n pseudorandom values of a random
%     variable x whose mean we want to estimate. fx is often defined as a
%     function of some random variable X with a simple distribution. The
%     input of fx should be a function handle, the output of fx should be
%     n function values.      
%
%     abstol --- the absolute error tolerance, which should be
%     positive, default value is 1e-2.
%
%     alpha --- the uncertainty, which should be a small positive
%     percentage. The default value is 1%.
%
%     nSig --- the number of samples used to compute the sample variance
%
%     fudge --- the standard deviation inflation factor
%
%                           Output Arguments
%
%     tmu --- the estimated value of the integral.
%
%     out_param.ntot --- total samples used.
%
%     out_param.var --- the sample variance.
%
%     out_param.time --- the time elapsed in seconds.
%

%This is a heuristic algorithm based on a Central Limit Theorem
%approximation combined with adaptive importance sampling.

if nargin < 5
   fudge = 1.2; %variance inflation factor
   if nargin < 4;
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 3
         alpha = 0.01; %uncertainty
         if nargin < 2
            abstol = 0.01; %absolute error tolerance
            if nargin < 1
                fx=@(x)max(x,0);%example function
            end
         end
      end
   end
end
nMax=1e8; %maximum number of samples allowed.
out_param.alpha = alpha; %save the input parameters to a structure
out_param.fudge = fudge;
out_param.nSig = nSig;

b=[-1,0,1];%estimated values for standard deviation.
tstart=tic; %starts the clock.
fx_b=@(t,b_value)fx(t+b_value).*exp(-t.*b_value-b_value.^2/2);
t_value=(randn(nSig,1));%generate the pseudorandom values drawn
% from the standard normal distribution

for i=1:numel(b)
    b_value=b(i);
    var_b(i)=var(fx_b(t_value,b_value));
end

[S_var,S_pos]=min(var_b);
out_param.var = S_var; %calculate the sample variance--stage 1
out_param.b_value=b(S_pos);%best variance
sig0 = sqrt(out_param.var); %standard deviation
sig0up = out_param.fudge.*sig0; %upper bound on the standard deviation
nmu = max(1,ceil((-norminv(alpha)*sig0up/abstol).^2)); 
   %number of samples needed for mean
assert(nmu<nMax,['nmu = ' int2str(nmu) ', which is too big']) 
   %don't exceed sample budget
tmu=mean(fx_b(randn(nmu,1),out_param.b_value)); %estimated mean
out_param.ntot=nSig+nmu; %total samples required
out_param.time=toc(tstart); %elapsed time    



end

