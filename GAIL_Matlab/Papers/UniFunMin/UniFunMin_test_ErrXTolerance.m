%This function describes the Experiment 3:  
% 
%  Generates Table 3.3 in the thesis with abstol=10^(-8), TolX=10^(-6), 
%  nrep=10000 and nmax=10^7
%
%  Xin Tong, A Guaranteed, Adaptive, Automatic Algorithm for Univatiate
%  Function Minimization, July 2014. 

%% Garbage collection and initialization
% clearvars -except testCase  %clear all variables except testCase
% close all %close all figures
function [tauvec,prob]=UniFunMin_test_ErrXTolerance(nrep,abstol,TolX,nmax)
tstart = tic;

%% Program parameters
in_param.abstol = abstol; %error tolerance
in_param.TolX = TolX;
in_param.nmax = nmax; %cost budget

%% Simulation parameters
n = nrep;
if (n >= 100)
    warning('off','MATLAB:funmin_g:exceedbudget');
    warning('off','MATLAB:funmin_g:peaky');
end;
a = 10.^(-4+3*rand(n,1));
z = 2.*a+(1-4*a).*rand(n,1);
tauvec = [11 101 1001]; %cone condition tau
ntau = length(tauvec);
ratio = 1./a;
exactmin = -1;
exactsolu = z;

%% Simulation
ntrapmat = zeros(nrep,ntau);
trueerrormat = ntrapmat;
truesolumat = ntrapmat;
newtaumat = ntrapmat;
tauchangemat = ntrapmat;
exceedmat = ntrapmat;
intnum = ntrapmat;

for i=1:ntau;
    for j=1:nrep;
        f = @(x) 0.5/a(j)^2*(-4*a(j)^2-(x-z(j)).^2-(x-z(j)-a(j)).*...
            abs(x-z(j)-a(j))+(x-z(j)+a(j)).*abs(x-z(j)+a(j))).*...
            (x>=z(j)-2*a(j)).*(x<=z(j)+2*a(j)); %test function
        in_param.nlo = (tauvec(i)+1)/2+1;
        in_param.nhi = in_param.nlo;
        [fmin,out_param] = funmin_g(f,in_param);
        ntrapmat(j,i) = out_param.npoints;
        newtaumat(j,i) = out_param.tau;
        estmin = fmin;
        trueerrormat(j,i) = abs(estmin-exactmin);
        tauchangemat(j,i) = out_param.tauchange;
        exceedmat(j,i) = out_param.exitflag;
        intnum(j,i) = size(out_param.intervals,2);
        for k=1:intnum(j,i)
            if exactsolu(j) <= out_param.intervals(2,k) && exactsolu(j) ... 
                    >= out_param.intervals(1,k) 
                truesolumat(j,i) = 1;
            end
        end
    end
end

prob.probinit = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); 
prob.probfinl = mean(repmat(ratio,1,ntau)<=newtaumat,1); 
prob.succnowarn=mean((trueerrormat<=in_param.abstol|truesolumat)&(~exceedmat),1); 
prob.succwarn=mean((trueerrormat<=in_param.abstol|truesolumat)&(exceedmat),1);    
prob.failnowarn=mean((trueerrormat>in_param.abstol&~truesolumat)&(~exceedmat),1);  
prob.failwarn=mean((trueerrormat>in_param.abstol&~truesolumat)&(exceedmat),1);  

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:ntau
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%% '],...
        [tauvec(i) 100*[prob.probinit(i) prob.probfinl(i) ...
        prob.succnowarn(i) prob.succwarn(i) prob.failnowarn(i)... 
        prob.failwarn(i)]])) 
end

%% Save output
gail.save_mat('UniFunMinOutput', 'ErrXToleranceTest',tauvec,prob,ntau);

toc(tstart)

warning('on','MATLAB:funmin_g:exceedbudget');
warning('on','MATLAB:funmin_g:peaky');

%% The following output was obtained on 2014-May
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.54%->21.00%   21.00%      0.00%   79.00%    0.00% 
%   101 33.42%->52.28%   52.28%      0.00%   47.72%    0.00% 
%  1001 66.15%->85.32%   85.33%      0.00%   14.67%    0.00% 
% Elapsed time is 1936.495643 seconds.

