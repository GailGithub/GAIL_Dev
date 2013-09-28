%Execution file of automatic guaranteed algorithm for function integration
%  Generates Table 3 in the paper
%
%  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell and Y. Zhang,
%  The Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, 
%  Not Balls, submitted for publication, arXiv.org:1303.2412 [math.NA]}, 
%  2013.

%% Preliminaries
format long e
clear all
close all
tstart = tic;

%% Program parameters
in_param.nmax=1e7; %maximum number of sample points
in_param.tol=1e-8; %error tolerance

%% Simulation parameters

% nrep = 100; %number of times to test, takes about a minute, can be changed 
nrep = 10000; %number of times to test used in the paper
if (nrep >= 1000)
    warning(' Need more than one hour to replicate the result in the paper! ')
    warning('off','MATLAB:integraltau_g:exceedbudget');
    warning('off','MATLAB:integraltau_g:peaky');
end;
avec = 10.^(rand(nrep,1).*3-4);
zvec = rand(nrep,1).*(1-4.*avec)+2.*avec;
x0vec = zvec-2*avec;
x1vec = zvec+2*avec;
tauvec = [10 100 1000]; %cone condition tau
ntau = length(tauvec);
ratio = 2./avec;

%% Simulation
% generate matrices to save data
Qmat=zeros(nrep,ntau);
npointsmat=Qmat;
errestmat=Qmat;
exactintmat=ones(nrep,ntau);
newtaumat=Qmat;
budgetmat=false(nrep,ntau);
tauchangemat=budgetmat;
timemat=Qmat;

% Qquadvec=zeros(nrep,1);
% Qintegralvec=Qquadvec;
% Qchebvec=Qquadvec;
% quadtimevec=Qquadvec;
% integraltimevec=Qquadvec;
% chebtimevec=Qquadvec;

% computes integrals for each function, each tau
for i=1:nrep
    if floor(i/100) == i/100
        display(i);
    end;
    %% Integrand
    a=avec(i);
    z=zvec(i);
    x0=x0vec(i);
    x1=x1vec(i);
    
    f = @(x) (1/(4*a^3)).*(4*a^2+(x-z).^2 ...
        +(x-z-a).*abs(x-z-a)...
        -(x-z+a).*abs(x-z+a)) ...
        .*(x>=x0).*(x<=x1); %test function

    %%Integrate
    for j=1:ntau
        in_param.tau=tauvec(j);
        tic
        [Q,out_param]=integraltau_g(f,in_param);
        timemat(i,j)=toc;
        Qmat(i,j)=Q;
        npointsmat(i,j)= out_param.npoints;
        newtaumat(i,j) = out_param.tau;
        budgetmat(i,j) = out_param.exceedbudget;
        errestmat(i,j)=out_param.errest;
        tauchangemat(i,j) = out_param.tauchange;
    end
    
    %computes integrals for each function using quad, integral, and
%     chebfun
%     tic;
%     Qquadvec(i)=quad(f,0,1,in_param.tol);
%     quadtimevec(i)=toc;
%     tic;
%     Qintegralvec(i)=integral(f,0,1,'AbsTol',in_param.tol);
%     integraltimevec(i)=toc;
%     tic;
%     Qchebvec(i)=sum(chebfun(f,[0,1]));
%     chebtimevec(i)=toc;
end

trueerrormat=abs(exactintmat-Qmat); % differences between true value and approximation of the integral

pini = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); %probability in initial cone
pfin = mean(repmat(ratio,1,ntau)<=newtaumat,1); %probability in final cone
exactsucc = mean(trueerrormat<=in_param.tol,1); %percentage of successful instants
succnowarn = mean((trueerrormat<=in_param.tol)&(~budgetmat),1); %percentage of successful instants for which the functions are in the cone
succwarn = mean((trueerrormat<=in_param.tol)&(budgetmat),1);    %percentage of successful instants for which the functions are not in the cone
failnowarn = mean((trueerrormat>in_param.tol)&(~budgetmat),1);  %percentage of failed instants for which the functions are in the cone
failwarn = mean((trueerrormat>in_param.tol)&(budgetmat),1); %percentage of failed instants for which the functions are not in the cone

% quadtrueerrorvec=abs(exactintmat(:,1)-Qquadvec);   % error of quad
% quadsuccessrate=mean(quadtrueerrorvec<=in_param.tol,1); % percentage of successful instants of quad
% quadavgtime=mean(quadtimevec,1);    % operation time of quad
% integraltrueerrorvec=abs(exactintmat(:,1)-Qintegralvec);   % error of integral
% integralsuccessrate=mean(integraltrueerrorvec<=in_param.tol,1); % percentage of successful instants of integral
% integralavgtime=mean(integraltimevec,1);    % operation time of integral
% chebtrueerrorvec=abs(exactintmat(:,1)-Qchebvec);   % error of chebfun
% chebsuccessrate=mean(chebtrueerrorvec<=in_param.tol,1); % percentage of successful instants of chebfun
% chebavgtime=mean(chebtimevec,1);    % operation time of chebfun


%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:ntau
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%% '],...
        [tauvec(i) 100*[pini(i) pfin(i) succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)]])) 
end
% display(['  success rate of quad = ', num2str(geomean(quadsuccessrate))])
% display(['  success rate of integral = ', num2str(geomean(integralsuccessrate))])
% display(['  success rate of chebfun = ', num2str(geomean(chebsuccessrate))])

%% Save Output
[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    'ConesPaperOutput',PATHNAMESEPARATOR','ConesPaperIntegralTest-',...
    datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.mat');
save(filename, ...
    'tauvec','pini','pfin','succnowarn', ...
    'succwarn','failnowarn','failwarn')
toc(tstart)
warning('on','MATLAB:integraltau_g:exceedbudget');
warning('on','MATLAB:integraltau_g:peaky');

%% The following output was obtained on 2013-August-05 by
%  from the data in
%       ConesPaperIntegrationTest02-Aug-2013-17-28-15.mat
%  by running the output section
%
%
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    10  0.00%->25.54%   25.35%      0.19%   74.46%    0.00% 
%   100 23.22%->58.04%   56.38%      1.66%   41.96%    0.00% 
%  1000 56.98%->87.61%   67.94%     19.84%   12.22%    0.00% 

