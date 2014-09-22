% This script is to test meanMCBernoulli_g algorithm
clear all;close all;clc
%format long e
in_param.errtype = 'rel';
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
in_param.alpha = 0.05;% default uncertainty
in_param.nmax = 1e10;
res = zeros(nrep,10);
for k = 1:nrep
    k
    p_true = 10^(-3+u(k)*2);
    in_param.abstol = 10^(-5+3*uu(k));
    in_param.reltol = 10^(-2+uu(k));
    Yrand=@(n) binornd(1,p_true,n,1);
    if strcmp(in_param.errtype,'abs')
        [p,out_param]=meanMCBernoulli_g(Yrand,in_param);
        % the results using meanMCBernoulli_g
        abserr = abs(p-p_true);
        %         numstr=horzcat((in_param.index), '       ',...
        %             num2str(in_param.abstol,'%10.5e'),'       ', num2str(abserr,'%10.5e'),'     ',...
        %             num2str(p,'%10.5e'), '       ', num2str(p_true,'%10.5e'),'      ', ...
        %             num2str(out_param.n,'%10.5e'),'      ', num2str(out_param.nmax));
        %         if abserr > in_param.abstol,% if error does not meet tolerance, mark it
        %             disp([numstr,'     ****']);
        %         else
        %             disp(numstr);
        %         end
        res(k,1) = in_param.abstol;
        res(k,2) = abserr;
        res(k,3) = p;
        res(k,4) = p_true;
        res(k,5) = out_param.n;
        res(k,6) = out_param.nmax;
        res(k,7) = abserr./in_param.abstol;
    end
    if strcmp(in_param.errtype,'rel')
        [p,out_param]=meanMCBernoulli_g(Yrand,in_param);% the results using cubMC_g
        abserr = abs(p-p_true);
        relerr = abs(p-p_true)/p;
        %         numstr=horzcat(num2str(in_param.index), '       ',...
        %             num2str(in_param.reltol),'       ', num2str(relerr,'%10.5e'),'     ',...
        %              num2str(p), '       ', num2str(p_true),'     ',...
        %             num2str(out_param.n),'      ', num2str(out_param.nmax));
        %         if relerr > in_param.reltol,% if error does not meet tolerance, mark it
        %             disp([numstr,'     ****']);
        %         else
        %             disp(numstr);
        %         end
        res(k,1) = in_param.reltol;
        res(k,2) = relerr;
        res(k,3) = p;
        res(k,4) = p_true;
        res(k,5) = out_param.n;
        res(k,6) = out_param.nmax;
        res(k,7) = relerr./in_param.reltol;
    end
%     if strcmp(in_param.index,'both')
%         [p,out_param]=meanMCBernoulli_g(Yrand,in_param);% the results using cubMC_g
%         relerr = abs(p-p_true)/p;
%         abserr = abs(p-p_true); 
%         res(k,1) = in_param.abstol;
%         res(k,2) = abserr;
%         res(k,3) = in_param.reltol;
%         res(k,4) = relerr;
%         res(k,5) = p;
%         res(k,6) = p_true;
%         res(k,7) = out_param.n;
%         res(k,8) = out_param.nmax;
%         res(k,9) = abserr./in_param.abstol;
%         res(k,10) = relerr./in_param.reltol;
% %         numstr=horzcat(num2str(in_param.index), '       ',...
% %             num2str(in_param.abstol),'       ', num2str(abserr),'     ',...
% %             num2str(in_param.reltol),'       ', num2str(relerr,'%10.5e'),'     ',...
% %             num2str(p), '       ', num2str(p_true),'     ',...
% %             num2str(out_param.n),'      ', num2str(out_param.nmax));
% %         if relerr > in_param.reltol || abserr > in_param.abstol;
% %             % if error does not meet tolerance, mark it
% %             disp([numstr,'     ****']);
% %         else
% %             disp(numstr);
% %         end
%     end   
%     res(1,:) = horzcat('abstol','abserr','p','p_true','n','nmax');
%     res(k,1) = in_param.abstol;
%     res(k,2) = abserr;
%     res(k,3) = p;
%     res(k,4) = p_true;
%     res(k,5) = out_param.n;
%     res(k,6) = out_param.nmax;
%     res(k,7) = abserr./in_param.abstol; 
end
timestamp=datestr(now);
timestamp(timestamp==' ')='_';
timestamp(timestamp==':')='.';
%end
%loglog(res(:,4),res(:,7),'r*')
filename = strcat('TestmeanMCBernoulli-on-', in_param.errtype ,'-',timestamp,'.mat');
save(filename)