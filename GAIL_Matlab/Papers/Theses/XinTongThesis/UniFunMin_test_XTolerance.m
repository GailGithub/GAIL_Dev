%This function describes the Experiment 2
% 
%  Generates Table 3.2 in the thesis abstol=0, TolX=10^(-6), nrep=10000 and
%  nmax=10^7
%
%  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for Univariate
%  Function Minimization. MS thesis, Illinois Institute of Technology,
%  2014.

function [tauvec,prob]=UniFunMin_test_XTolerance(nrep,TolX,nmax)
tstart = tic;

%% Program parameters
in_param.abstol = 0; %error tolerance
in_param.TolX = TolX;
in_param.nmax = nmax; %cost budget

%% Simulation parameters
n = nrep;
if (n >= 100)
    warning('off','GAIL:funmin01_g:exceedbudget');
    warning('off','GAIL:funmin01_g:peaky');
end;
a = 10.^(-4+3*rand(n,1));
z = 2.*a+(1-4*a).*rand(n,1);
tauvec = [11 101 1001]; %cone condition tau
ntau = length(tauvec);
ratio = 1./a;
exactsolu = z;

%% Simulation
ntrapmat = zeros(nrep,ntau);
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
        in_param.ninit = (tauvec(i)+1)/2+1;
        [fmin,out_param] = funmin01_g(f,in_param);
        ntrapmat(j,i) = out_param.npoints;
        newtaumat(j,i) = out_param.tau;
        tauchangemat(j,i) = out_param.tauchange;
        exceedmat(j,i) = out_param.exceedbudget;
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
prob.succnowarn = mean((truesolumat)&(~exceedmat),1); 
prob.succwarn = mean((truesolumat)&(exceedmat),1);    
prob.failnowarn = mean((~truesolumat)&(~exceedmat),1);   
prob.failwarn = mean((~truesolumat)&(exceedmat),1);

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
gail.save_mat('UniFunMinOutput', 'XToleranceTest',true,tauvec,prob,ntau)

toc(tstart)

warning('on','GAIL:funmin01_g:exceedbudget');
warning('on','GAIL:funmin01_g:peaky');

%% The following output was obtained on 2015-Feb
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.28%->20.56%   20.56%      0.00%   79.44%    0.00% 
%   101 33.02%->52.20%   52.20%      0.00%   47.80%    0.00% 
%  1001 65.90%->84.72%   80.62%      4.10%   15.28%    0.00% 

