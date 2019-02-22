%CONEPAPER_TEST_FUNAPPX_G Generate Table 3. in Cones not ball paper Run automatic guaranteed algorithm for function approximation
%  Generates Table 3 in the paper
%
%  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell and Y. Zhang,
%  The Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, 
%  Not Balls, submitted for publication, arXiv.org:1303.2412 [math.NA]}, 
%  2013.
%

%% Preliminaries
%clear all, close all
%clearvars -except testCase
function [succnowarn,succwarn]=conepaper_test_funappx_g(nrep,nmax,abstol)
tstart = tic;

%% Program parameters
in_param.abstol = abstol; %error tolerance
in_param.nmax = nmax; %cost budget

%% Simulation parameters
%nrep = 100; %number of times to test, takes about a minute, can change
%nrep = 10000; %number of times to test used in the paper
if (nrep >= 1000)
    warning(' Need more than one hour to replicate the result in the paper! ')
    warning('off','GAIL:funappxtau_g:exceedbudget');
    warning('off','GAIL:funappxtau_g:peaky');
end;
a = 10.^(rand(nrep,1).*3-4);
z = rand(nrep,1).*(1-4*a)+2*a;
x0 = z-2*a;
x1 = z+2*a;
tauvec = [10 100 1000]; %cone condition tau
ntau = length(tauvec);
ratio = 1./a;
gnorm = 1./a;

%% Simulation
ntrapmat = zeros(nrep,ntau);
trueerrormat = ntrapmat;
newtaumat = ntrapmat;
timemat = ntrapmat;
tauchangemat = ntrapmat;
exceedmat = ntrapmat;
ballmat = ntrapmat;

for i=1:ntau;
    for j=1:nrep;
        f = @(x) 1/(2*a(j)^2)*(4*a(j)^2+(x-z(j)).^2+...
        (x-z(j)-a(j)).*abs(x-z(j)-a(j))-(x-z(j)+a(j)).*abs(x-z(j)+a(j))).*...
        (x>=x0(j)).*(x<=x1(j)); %test function
        in_param.tau = tauvec(i);
        tic
        [fappx,out_param] = funappxtau_g(f,in_param);
        chebfun
        timemat(j,i) = toc;
        ntrapmat(j,i) = out_param.npoints;
        newtaumat(j,i) = out_param.tau;
        xx = [rand(1,50000) rand(1,50000)*(4*a(j))+(z(j)-2*a(j))];
        yy1 = fappx(xx);
        exactyy = f(xx);
        trueerrormat(j,i) = max(abs(yy1-exactyy));
        tauchangemat(j,i) = out_param.tauchange;
        exceedmat(j,i) = out_param.exceedbudget;
        ballmat(j,i) = out_param.ballradius;
    end;
%     display(['When tau = ', num2str(tau(i))]) 
%     display(['  probability of f in the initial cone = ', num2str(pini)])
%     display(['  probability of f in the final cone = ', num2str(pfin)])
%     display(['  success rate of algorithm = ', num2str(exactsucc)])
%     display(['  success, no warning = ', num2str(succnowarn)])
%     display(['  success, warning = ', num2str(succwarn)])
%     display(['  failure, no warning = ', num2str(failnowarn)])
%     display(['  failure, warning = ', num2str(failwarn)])
%     display(['  average time cost = ', num2str(geomean(timemat(:,i)))])
end;

pini = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); %probability in initial cone
pfin = mean(repmat(ratio,1,ntau)<=newtaumat,1); %probability in final cone
exactsucc = mean(trueerrormat<=in_param.abstol,1); %percentage of successful instants
succnowarn = mean((trueerrormat<=in_param.abstol)&(~exceedmat),1); %percentage of successful instants for which the functions are in the cone
succwarn = mean((trueerrormat<=in_param.abstol)&(exceedmat),1);    %percentage of successful instants for which the functions are not in the cone
failnowarn = mean((trueerrormat>in_param.abstol)&(~exceedmat),1);  %percentage of failed instants for which the functions are in the cone
failwarn = mean((trueerrormat>in_param.abstol)&(exceedmat),1);

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:3
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%% '],...
        [tauvec(i) 100*[pini(i) pfin(i) succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)]])) 
end

%% Save Output
% [GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
% filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
%     'ConesPaperOutput',PATHNAMESEPARATOR','ConesPaperFunAppxTest-',...
%     datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.mat');
% save(filename, ...
%     'tauvec','pini','pfin','succnowarn', ...
%     'succwarn','failnowarn','failwarn')
gail.save_mat('ConesPaperOutput', 'ConesPaperFunAppxTest', true, tauvec,pini,...
    pfin,succnowarn,succwarn,failnowarn,failwarn);
toc(tstart)
warning('on','GAIL:funappxtau_g:exceedbudget');
warning('on','GAIL:funappxtau_g:peaky');
end
%% The following output was obtained on 2013-Sep-03 by
%  from the data in
%       ConesPaperFunAppxTest-03-Sep-2013-18-37-05.mat
%  by running the output section
%
%
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%   10  0.00%->26.12%   25.75%      0.24%   73.88%    0.13% 
%  100 32.95%->57.74%   56.30%      0.76%   42.26%    0.68% 
% 1000 66.59%->88.21%   75.82%      4.68%   11.78%    7.72% 
%% The following output was obtained on 2015-Jan-24 by
%  from the data in
%       ConesPaperFunAppxTest-2015-Jan-24-13-34-24.mat
%  by running the output section
%
%
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    10  0.00%->25.38%   25.13%      0.17%   74.62%    0.08% 
%   100 32.84%->57.65%   56.36%      0.62%   42.35%    0.67% 
%  1000 66.70%->88.62%   75.71%      4.52%   11.38%    8.39% 
