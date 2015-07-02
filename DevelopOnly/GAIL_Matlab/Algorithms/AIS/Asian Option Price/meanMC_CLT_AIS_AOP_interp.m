function [tmu,out_param]=meanMC_CLT_AIS_AOP_interp(d,So,K, r, sigma, T,abstol,alpha,nSig,fudge)
%MEANMC_CLT Monte Carlo method to estimate the mean of a random variable
%
%   Estimates the mean, of a random fuction to
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
%     d --- number of dimensions.
%
%     So --- Initial stock price.
%     
%     K ---
%
%     r ---
%
%     sigma ---
%
%     T --- time.
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

if nargin < 10
   fudge = 1.2; %variance inflation factor
   if nargin < 9;
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 8
         alpha = 0.01; %uncertainty
         if nargin < 7
            abstol = 0.002; %absolute error tolerance
            if nargin < 6
                T = 1;
                if nargin < 5
                    sigma = 0.5;
                    if nargin < 4
                        r = 0;
                        if nargin < 3
                            K = 15;
                            if nargin < 2
                                So = 10;
                                if nargin < 1
                                    d = 2;
                                end
                            end
                        end
                    end
                end 
            end
         end
      end
   end
end
nMax=1e8; %maximum number of samples allowed.
out_param.alpha = alpha; %save the input parameters to a structure
out_param.fudge = fudge;
out_param.nSig = nSig;

tstart=tic; %starts the clock.

b=[-4,0,0.5]./d;%estimated values for standard deviation.

%z is a d dimensional normal random vector
%b is the shift in the mean
deltaT = T./d;

SK=@(z,bval) max(mean(So.*exp(cumsum((r-((sigma.^2)./2)).*deltaT + ...
   sigma.*sqrt(deltaT).*(z-bval),2)),2)-K,0);
% 
% for j=1:d; 
%     tj=j.*T./d;
%     %sum of S-K:
%     SK=@(z,b)(So.*exp(r-((sigma.^2)./2)).*tj + sigma.*sqrt(T./d).*z-b)-K;
% end

%integrand:
fx=@(z,bval) SK(z,bval).*exp(-r.*T).*exp(sum(z,2).*bval-(d.*(bval.^2)./2)); 

z=randn(nSig,d);%generate the pseudorandom values drawn
% from the standard normal distribution
%computes each value of tj
var_b=b;
for i=1:numel(b)
    b_value=b(i);    
    var_b(i)=var(fx(z,b_value));
end
var_b

%parabolic interpolation:
A=[b(1).^2 b(1) 1; b(2).^2 b(2) 1; b(3).^2 b(3) 1];
p=A\var_b';
fmin=@(x)p(3)+p(2)*x+p(1)*(x.^2);
[x]=fminbnd(fmin,b(1),b(3))
var_bx=var(fx(z,x))

%%
% 
% <html>
% <table border=1><tr><td>one</td><td>two</td></tr></table>
% </html>
% 
[S_var,S_pos]=min(var_b);
if var_bx < S_var
    out_param.b_value = x;
    out_param.var = var_bx;
else
    out_param.b_value = b(S_pos);
    out_param.var = S_var;
end
%out_param.var = min(S_var, var_bx); %calculate the sample variance--stage 1
%out_param.b_value = x;%best variance
sig0 = sqrt(out_param.var); %standard deviation
sig0up = out_param.fudge.*sig0; %upper bound on the standard deviation
nmu = max(1,ceil((-norminv(alpha)*sig0up/abstol).^2)); 
   %number of samples needed for mean
assert(nmu<nMax,['nmu = ' int2str(nmu) ', which is too big']) 
   %don't exceed sample budget
tmu=mean(fx(randn(nmu,1),out_param.b_value)); %estimated mean
out_param.ntot=nSig+nmu; %total samples required
out_param.time=toc(tstart); %elapsed time    



end

