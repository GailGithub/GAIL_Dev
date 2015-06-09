%This function describes the Experiment 1
% 
%  Generates Table 3.1 in the thesis with abstol=10^(-8), TolX=0, 
%  and nrep=10000 nmax=10^7
%
%  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for Univariate
%  Function Minimization. MS thesis, Illinois Institute of Technology,
%  2014.

function [tauvec,prob]=UniFunMin_test_ErrTolerance(nrep,abstol,nmax)
tstart = tic;

%% Program parameters
in_param.abstol = abstol; %error tolerance
in_param.TolX = 0;
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
exactmin = -1;

%% Simulation
ntrapmat = zeros(nrep,ntau);
trueerrormat = ntrapmat;
newtaumat = ntrapmat;
tauchangemat = ntrapmat;
exceedmat = ntrapmat;

for i=1:ntau;
    for j=1:nrep;
        f = @(x) 0.5/a(j)^2*(-4*a(j)^2-(x-z(j)).^2-(x-z(j)-a(j)).*...
            abs(x-z(j)-a(j))+(x-z(j)+a(j)).*abs(x-z(j)+a(j))).*...
            (x>=z(j)-2*a(j)).*(x<=z(j)+2*a(j)); %test function
        in_param.ninit = (tauvec(i)+1)/2+1;
        [fmin,out_param] = funmin01_g(f,in_param);
        ntrapmat(j,i) = out_param.npoints;
        newtaumat(j,i) = out_param.tau;
        estmin = fmin;
        trueerrormat(j,i) = abs(estmin-exactmin);
        tauchangemat(j,i) = out_param.tauchange;
        exceedmat(j,i) = out_param.exceedbudget;
    end
end

prob.probinit = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); 
prob.probfinl = mean(repmat(ratio,1,ntau)<=newtaumat,1); 
prob.succnowarn=mean((trueerrormat<=in_param.abstol)&(~exceedmat),1); 
prob.succwarn=mean((trueerrormat<=in_param.abstol)&(exceedmat),1);    
prob.failnowarn=mean((trueerrormat>in_param.abstol)&(~exceedmat),1);  
prob.failwarn=mean((trueerrormat>in_param.abstol)&(exceedmat),1);  

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
gail.save_mat('UniFunMinOutput', 'ErrToleranceTest',true,tauvec,prob,ntau);

toc(tstart)

warning('on','GAIL:funmin01_g:exceedbudget');
warning('on','GAIL:funmin01_g:peaky');

%% The following output was obtained on 2015-Feb
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.30%->56.85%   47.62%      3.58%   43.13%    5.67% 
%   101 33.39%->75.11%   66.35%      3.42%   24.88%    5.35% 
%  1001 66.69%->93.29%   81.37%      5.72%    6.70%    6.21% 
