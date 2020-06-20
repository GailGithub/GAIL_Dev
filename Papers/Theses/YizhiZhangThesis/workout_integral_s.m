%WORKOUT_INTEGRAL_G Calls automatic guaranteed algorithm for univariate integration

function [succnowarn, succwarn]=workout_integral_s(nrep,nmax,abstol)
%% Garbage cleanup
format long e
close all;
tstart = tic;
warning('off','GAIL:integral_g:exceedbudget');
warning('off','GAIL:integral_g:spiky');
warning off
%% Program parameters
% in_param.nmax=1e7; %maximum number of sample points
% in_param.abstol=1e-8; %error tolerance

%% Simulation parameters

%nrep = 100; %number of times to test, takes about a minute 
deltavec = 10.^(rand(nrep,1).*3-4);
tvec = rand(nrep,1).*(1-4.*deltavec);
%x0vec = tvec-2*deltavec;
%delta4vec = tvec+4*deltavec;

hcuttrapvec = [0.1 0.01 0.001];
% hcutsimpvec = hcuttrapvec*3;
hcutsimpvec = hcuttrapvec;
hcutvec=[hcuttrapvec hcutsimpvec];
ninittrapvec = ceil(2./hcuttrapvec)+2;
ninitsimpvec = ceil(6./hcutsimpvec)+2;
ninitvec = [ninittrapvec ninitsimpvec]; %starting number of points
nninit = length(ninitvec);


% avec = 10.^(rand(nrep,1).*3-4);
% zvec = rand(nrep,1).*(1-4.*avec)+2.*avec;
% x0vec = zvec-2*avec;
% x1vec = zvec+2*avec;
% ninitvec = [22 202 2002]; %starting number of points
% hcutvec = [0.1 0.01 0.001];
% nninit = length(ninitvec);
% ratio = 2./avec;

%% Simulation
% generate matrices to save data
Qmat=zeros(nrep,nninit);
npointsmat=Qmat;
errestmat=Qmat;
exactintmat=ones(nrep,nninit);
newhcutmat=Qmat;
budgetmat=false(nrep,nninit);
conechangemat=budgetmat;
timemat=Qmat;
integralmat=zeros(nrep,1);
chebfunmat=integralmat;

% computes integrals for each function, each ninit
for i=1:nrep
    if floor(i/100) == i/100
        disp(i);
    end
    %% Integrand
    t=tvec(i);
    delta=deltavec(i)/4;
    
    f = @(x) (1/delta.^4).*((x-t).^3/6.*(x>=t).*(x<t+delta)...
        +(-3.*(x-t).^3+12.*delta.*(x-t).^2-12.*delta.^2.*(x-t)+4.*delta.^3)/6.*(x>=t+delta).*(x<t+2*delta)...
        +(3.*(x-t).^3-24.*delta.*(x-t).^2+60.*delta.^2.*(x-t)-44.*delta.^3)/6.*(x>=t+2*delta).*(x<t+3*delta)...
        +(t+4.*delta-x).^3/6.*(x>=t+3*delta).*(x<=t+4*delta)); %test function

    %%Integrate
    for j=1:nninit/2
        hcut=hcutvec(j);
        tic
        [q,out_param]=integral_t(f,'nmax',nmax,'abstol',abstol,'a',0,'b',1,'hcut',hcut);
        timemat(i,j)=toc;
        Qmat(i,j)=q;
        npointsmat(i,j)= out_param.npoints;
        newhcutmat(i,j) = out_param.hcut;
        budgetmat(i,j) = out_param.exceedbudget;
        errestmat(i,j) = out_param.errest;
        conechangemat(i,j) = out_param.conechange;
    end
    
    for j=nninit/2+1:nninit
        hcut=hcutvec(j);
        tic
        [q,out_param]=integral_s(f,'nmax',nmax,'abstol',abstol,'a',0,'b',1,'hcut',hcut);
        timemat(i,j)=toc;
        Qmat(i,j)=q;
        npointsmat(i,j)= out_param.npoints;
        newhcutmat(i,j) = out_param.hcut;
        budgetmat(i,j) = out_param.exceedbudget;
        errestmat(i,j) = out_param.errest;
        conechangemat(i,j) = out_param.conechange;
    end
    integralmat(i)=integral(f,0,1);
    chebfunmat(i)=sum(chebfun(f,[0,1]));
end

trueerrormat=abs(exactintmat-Qmat); % differences between true value and approximation of the integral
trueerrorintegralmat=abs(exactintmat(:,1)-integralmat);
trueerrorchebfunmat=abs(exactintmat(:,1)-chebfunmat);

% exactsucc = mean(trueerrormat<=abstol,1); %percentage of successful instants
succnowarn = mean((trueerrormat<=abstol)&(~conechangemat),1); %percentage of successful instants for which the functions are in the cone
succwarn = mean((trueerrormat<=abstol)&(conechangemat),1);    %percentage of successful instants for which the functions are not in the cone
failnowarn = mean((trueerrormat>abstol)&(~conechangemat),1);  %percentage of failed instants for which the functions are in the cone
failwarn = mean((trueerrormat>abstol)&(conechangemat),1); %percentage of failed instants for which the functions are not in the cone
succintegral = mean(trueerrorintegralmat<=abstol);
succchebfun = mean(trueerrorchebfunmat<=abstol);

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
disp(' ')
disp('          Success    Success   Failure    Failure')
disp(' ninit   No Warning  Warning  No Warning  Warning')
for i=1:nninit
    fprintf(['%5.0f %10.2f%%' ...
        '%10.2f%% %8.2f%% %8.2f%% \n'],...
        [ninitvec(i) 100*[succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)]])
end
disp(['integral= ' num2str(100*succintegral)])
disp(['chebfun= ' num2str(100*succchebfun)]) 

%% Save Output

time = toc(tstart);

gail.save_mat('WorkoutIntegralOutput', 'WorkoutIntegralTest',true, nrep,time,...
        succnowarn,succwarn,failnowarn,failwarn);
warning on
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