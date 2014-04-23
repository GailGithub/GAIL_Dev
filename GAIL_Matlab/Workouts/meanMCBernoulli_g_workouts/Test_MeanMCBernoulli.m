% This script is to test meanMCBernoulli_g algorithm
clear all;close all;clc
%format long e
in_param.index = 'abs-clt';
% if strcmp(in_param.index,'abs')
%     disp(horzcat('index       abstol            abserr          p_hat             p_true             n           nmax '));
%     disp(        '-------------------------------------------------------------------------------------------------------');
% end
% if strcmp(in_param.index,'rel')
%     disp(horzcat('index      reltol     relerr       p_hat         p_true       n      nmax'));
%     disp(        '----------------------------------------------------------------------');
% end
% if strcmp(in_param.index,'both')
%     disp(horzcat('index      abstol     abserr     reltol     relerr      p_hat         p_true       n      nmax'));
%     disp(        '----------------------------------------------------------------------');
% end
nrep = 500;
u = rand(nrep,1);
uu = rand(nrep,1);
in_param.alpha = 1e-2;% default uncertainty
in_param.nmax = 1e9;
res = zeros(nrep,10);
for k = 1:nrep
    p_true = 10^(-3+u(k)*2);
    in_param.abstol = p_true*10^(-2+uu(k));
    %in_param.reltol = 10^(-2+uu(k));
    Yrand=@(n) binornd(1,p_true,n,1);

        [p,out_param]=meanMCBernoulli(Yrand,in_param);% the results using cubMC_g
        abserr = abs(p-p_true);
        res(k,1) = in_param.abstol;
        res(k,2) = abserr;
        res(k,3) = p;
        res(k,4) = p_true;
        res(k,5) = out_param.n_clt;
        res(k,6) = out_param.nmax;
        res(k,7) = abserr./in_param.abstol;
    end
%     if strcmp(in_param.index,'rel')
%         [p,out_param]=meanMCBernoulli_g(Yrand,in_param);% the results using cubMC_g
%         abserr = abs(p-p_true);
%         relerr = abs(p-p_true)/p;
%         res(k,1) = in_param.reltol;
%         res(k,2) = relerr;
%         res(k,3) = p;
%         res(k,4) = p_true;
%         res(k,5) = out_param.n;
%         res(k,6) = out_param.nmax;
%         res(k,7) = relerr./in_param.reltol;
%     end
%     if strcmp(in_param.index,'both')
%         [p,out_param]=meanMCBernoulli_g(Yrand,in_param);% the results using cubMC_g
%         relerr = abs(p-p_true)/p;
%         abserr = abs(p-p_true); 
%         res(k,1) = in_param.reltol;
%         res(k,2) = relerr;
%         res(k,3) = in_param.reltol;
%         res(k,4) = relerr;
%         res(k,5) = p;
%         res(k,6) = p_true;
%         res(k,7) = out_param.n;
%         res(k,8) = out_param.nmax;
%         res(k,9) = abserr./in_param.abstol;
%         res(k,10) = relerr./in_param.reltol;
%     end
timestamp=datestr(now);
timestamp(timestamp==' ')='_';
timestamp(timestamp==':')='.';
%end
%loglog(res(:,4),res(:,7),'r*')
filename = strcat('TestmeanMCBernoulli-on-', in_param.index ,'-',timestamp,'.mat');
save(filename)