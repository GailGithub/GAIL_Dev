function [tmu,out_param]=meanMC_CLT_AIS_AOP_g(Y,b_vec,abstol,alpha,nSig,fudge)
%   MEANMC_CLT_AIS_AOP_PutCall uses adaptive importance sampling 
%   to estimates the asian option price within a specified error tolerance, 
%   i.e., | mu - tmu | <= abstol with probability at least 1-alpha, where 
%   abstol is the absolute error tolerance using the adaptive importance sampling.  
%
%
%                           Input Arguments
%
%
%     B --- Vector with three values for standard deviation.
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

if nargin < 6
   fudge = 1.2; %variance inflation factor
   if nargin < 5
      nSig = 1e4; %number of samples to estimate variance
      if nargin < 4
         alpha = 0.01; %uncertainty
         if nargin < 3
            abstol = 0.002; %absolute error tolerance
            if nargin < 2
                b_vec=[-1,0,1];   
                if nargin < 1
                   Y1=@(x,b)max(x-b,0).*exp(x.*b-(b.^2)./2);
                   Y=@(n,b)Y1(randn(n,1),b);
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

%c_p=input('\n Type "1" for PUT or "2" for CALL:\n\n');

tstart=tic; %starts the clock.


var_b=b_vec;
for i=1:numel(b_vec)
    b_value=b_vec(i);    
    var_b(i)=var(Y(nSig,b_value));
end
var_b

%parabolic interpolation:
A=[b_vec'.^2 b_vec' ones(3,1)];
p=A\var_b';
fmin=@(x)p(3)+p(2)*x+p(1)*(x.^2);
[x]=fminbnd(fmin,b_vec(1),b_vec(3));
var_bx=var(Y(nSig,x));
var_bx
x


[S_var,S_pos]=min(var_b);
if var_bx < S_var
    out_param.b_value = x;
    out_param.var = var_bx;
else
    out_param.b_value = b_vec(S_pos);
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
tmu=mean(Y(nmu,out_param.b_value)); %estimated mean
out_param.ntot=nSig+nmu; %total samples required
out_param.time=toc(tstart); %elapsed time    



end

