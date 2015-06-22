% This script is to test meanMCBer_g algorithm
%clear all;close all;clc
function [ut_abserr, ut_abstol, nrep, res, out_param]=Test_meanMCBer_g
format long e
disp(horzcat(' abstol             abserr            pHat                p               n               nmax '));
disp(        '-------------------------------------------------------------------------------------------------------');
%nrep = 500;%in paper, we use 500 replications
nrep = 50;
u = rand(nrep,1);
uu = rand(nrep,1);
in_param.alpha = 0.05;% default uncertainty
in_param.nmax = 1e10;
res = zeros(nrep,10);
ut_abserr = zeros(nrep,1);
ut_abstol = zeros(nrep,1);
warning('off','GAIL:meanMCBer_g:nabsexceednmax')
for k = 1:nrep
    p = 10^(-3+u(k)*2);
    in_param.abstol = 10^(-5+3*uu(k));
    Yrand=@(n) binornd(1,p,n,1);
    [pHat,out_param]=meanMCBer_g(Yrand,in_param);
    % the results using meanMCBernoulli_g
    abserr = abs(pHat-p);
    ut_abserr(k) = abserr;
    ut_abstol(k) = in_param.abstol;
    numstr=horzcat(num2str(in_param.abstol,'%10.5e'),'       ', num2str(abserr,'%10.5e'),'     ',...
        num2str(pHat,'%10.5e'), '       ', num2str(p,'%10.5e'),'      ', ...
        num2str(out_param.n,'%10.5e'),'      ', num2str(out_param.nmax));
    if k <= 10
        if abserr > in_param.abstol,% if error does not meet tolerance, mark it
            disp([numstr,'     ****']);
        else
            disp(numstr);
        end
    end
    res(k,1) = in_param.abstol;
    res(k,2) = abserr;
    res(k,3) = pHat;
    res(k,4) = p;
    res(k,5) = out_param.n;
    res(k,6) = out_param.nmax;
    res(k,7) = abserr./in_param.abstol;
end
gail.save_mat('meanMCBerPaperOutput','TestmeanMCBer-on-abs-',true,ut_abstol,ut_abserr,res,nrep,out_param)
warning('on','GAIL:meanMCBer_g:nabsexceednmax')
end
