%TRYOUT_INTEGRAL_G_FJH  Try out adapttrap

%put the paper name here including the name of the figure or table that it
%produces

%% Garbage cleanup
format long e
clear all
close all
tstart = tic;

%% Program parameters
in_param.nmax=1e7; %maximum number of sample points
in_param.tol=1e-8; %error tolerance

%% Simulation parameters

nrep = 100; %number of times to test, takes about a minute, can be changed 
%nrep = 10000; %number of times to test used in the paper
avec = 10.^(rand(nrep,1).*3-4);
zvec = rand(nrep,1).*(1-4.*avec)+2.*avec;
x0vec = zvec-2*avec;
x1vec = zvec+2*avec;
tauvec = [10 100 1000]; %cone condition tau
ntau = length(tauvec);
ratio = 2./avec;

%% Simulation
Qmat=zeros(nrep,ntau);
npointsmat=Qmat;
errestmat=Qmat;
exactintmat=ones(nrep,ntau);
newtaumat=Qmat;
budgetmat=false(nrep,ntau);
tauchangemat=budgetmat;
timemat=Qmat;

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
end

%% Summarize Output
trueerrormat=abs(exactintmat-Qmat);
pini = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); %probability in initial cone
pfin = mean(repmat(ratio,1,ntau)<=newtaumat,1); %probability in final cone
exactsucc = mean(trueerrormat<=in_param.tol,1);
succnowarn = mean((trueerrormat<=in_param.tol)&(~budgetmat),1);
succwarn = mean((trueerrormat<=in_param.tol)&(budgetmat),1);
failnowarn = mean((trueerrormat>in_param.tol)&(~budgetmat),1);
failwarn = mean((trueerrormat>in_param.tol)&(budgetmat),1);

%% Save Output
save(['ConesPaperIntegrationTest' datestr(now,'dd-mmm-yyyy-HH-MM-SS') '.mat'])

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:ntau
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%%'],...
        [tauvec(i) 100*[pini(i) pfin(i) succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)]])) 
end

toc(tstart)

%% The following output was obtained on 2013-August-05 by
%  from the data in
%       ConesPaperIntegrationTest02-Aug-2013-17-28-15.mat
%  by running the output section
%  Difference of 1% are to be expected if n=10000
%
%
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    10  0.00%->25.54%   25.35%      0.19%   74.46%    0.00% 
%   100 23.22%->58.04%   56.38%      1.66%   41.96%    0.00% 
%  1000 56.98%->87.61%   67.94%     19.84%   12.22%    0.00% 
