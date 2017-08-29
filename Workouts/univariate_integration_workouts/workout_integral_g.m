%WORKOUT_INTEGRAL_G Calls automatic guaranteed algorithm for univariate integration

function [succnowarn, succwarn,pfin]=workout_integral_g(nrep,nmax,abstol)
%% Garbage cleanup
format long e
close all;
tstart = tic;
warning('off','GAIL:integral_g:exceedbudget');
warning('off','GAIL:integral_g:peaky');

%% Program parameters
% in_param.nmax=1e7; %maximum number of sample points
% in_param.abstol=1e-8; %error tolerance

%% Simulation parameters

%nrep = 100; %number of times to test, takes about a minute 
avec = 10.^(rand(nrep,1).*3-4);
zvec = rand(nrep,1).*(1-4.*avec)+2.*avec;
x0vec = zvec-2*avec;
x1vec = zvec+2*avec;
ninitvec = [7 52 502]; %starting number of points
tauvec = [10 100 1000];
nninit = length(ninitvec);
ratio = 2./avec;

%% Simulation
% generate matrices to save data
Qmat=zeros(nrep,nninit);
npointsmat=Qmat;
errestmat=Qmat;
exactintmat=ones(nrep,nninit);
newtaumat=Qmat;
budgetmat=false(nrep,nninit);
tauchangemat=budgetmat;
timemat=Qmat;

% computes integrals for each function, each ninit
for i=1:nrep
%     if floor(i/100) == i/100
%         display(i);
%     end;
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
    for j=1:nninit
        ninit=ninitvec(j);
        tic
        [q,out_param]=integral_g(f,'nmax',nmax,'abstol',abstol,'a',0,'b',1,'nlo',ninit,'nhi',ninit);
        timemat(i,j)=toc;
        Qmat(i,j)=q;
        npointsmat(i,j)= out_param.npoints;
        newtaumat(i,j) = out_param.tau;
        budgetmat(i,j) = out_param.exceedbudget;
        errestmat(i,j) = out_param.errest;
        %tauchangemat(i,j) = out_param.tauchange;
    end
end

trueerrormat=abs(exactintmat-Qmat); % differences between true value and approximation of the integral

pini = mean(repmat(ratio,1,nninit)<=repmat(tauvec,nrep,1),1); %probability in initial cone
pfin = mean(repmat(ratio,1,nninit)<=newtaumat,1); %probability in final cone
% exactsucc = mean(trueerrormat<=abstol,1); %percentage of successful instants
succnowarn = mean((trueerrormat<=abstol)&(~budgetmat),1); %percentage of successful instants for which the functions are in the cone
succwarn = mean((trueerrormat<=abstol)&(budgetmat),1);    %percentage of successful instants for which the functions are not in the cone
failnowarn = mean((trueerrormat>abstol)&(~budgetmat),1);  %percentage of failed instants for which the functions are in the cone
failwarn = mean((trueerrormat>abstol)&(budgetmat),1); %percentage of failed instants for which the functions are not in the cone

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' ninit    In Cone    No Warning  Warning No Warning Warning')
for i=1:nninit
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%% '],...
        [ninitvec(i) 100*[pini(i) pfin(i) succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)]])) 
end
 
%% Save Output

time = toc(tstart);

gail.save_mat('WorkoutIntegralOutput', 'WorkoutIntegralTest',true, nrep,time,...
        succnowarn,succwarn,failnowarn,failwarn);

warning('on','GAIL:integral_g:exceedbudget');
warning('on','GAIL:integral_g:peaky');

%% The following output was obtained on 2013-August-06 by
%  from the data in
%       ConesPaperIntegrationTest02-Aug-2013-11-15-16.mat
%  by running the output section
%
%
%         Probability    Success   Success   Failure  Failure
%  ninit    In Cone    No Warning  Warning No Warning Warning
%    10  0.00%->25.37%   25.18%      0.19%   74.63%    0.00% 
%   100 23.01%->57.56%   55.81%      1.75%   42.44%    0.00% 
%  1000 56.81%->88.45%   68.06%     20.54%   11.40%    0.00% 

end