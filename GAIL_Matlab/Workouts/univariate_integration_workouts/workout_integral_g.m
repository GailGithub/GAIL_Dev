%WORKOUT_INTEGRAL_G Calls automatic guaranteed algorithm for univariate integration

function [succnowarn, succwarn]=workout_integral_g(nrep,nmax,abstol)
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
hcutvec = [0.1 0.01 0.001];
ninitvec = [62 602 6002]; %starting number of points
nninit = length(ninitvec);



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



% computes integrals for each function, each ninit
for i=1:nrep
%     if floor(i/100) == i/100
%         disp(i);
%     end
    %% Integrand
    t=tvec(i);
    delta=deltavec(i)/4;
    
    f = @(x) (1/delta.^4).*((x-t).^3/6.*(x>=t).*(x<t+delta)...
        +(-3.*(x-t).^3+12.*delta.*(x-t).^2-12.*delta.^2.*(x-t)+4.*delta.^3)/6.*(x>=t+delta).*(x<t+2*delta)...
        +(3.*(x-t).^3-24.*delta.*(x-t).^2+60.*delta.^2.*(x-t)-44.*delta.^3)/6.*(x>=t+2*delta).*(x<t+3*delta)...
        +(t+4.*delta-x).^3/6.*(x>=t+3*delta).*(x<=t+4*delta)); %test function

    %%Integrate
    for j=1:nninit
        ninit=ninitvec(j);
        tic
        [q,out_param]=integral_g(f,'nmax',nmax,'abstol',abstol,'a',0,'b',1,'nlo',ninit,'nhi',ninit);
        timemat(i,j)=toc;
        Qmat(i,j)=q;
        npointsmat(i,j)= out_param.npoints;
        newhcutmat(i,j) = out_param.hcut;
        budgetmat(i,j) = out_param.exceedbudget;
        errestmat(i,j) = out_param.errest;
        conechangemat(i,j) = out_param.conechange;
    end
end

trueerrormat=abs(exactintmat-Qmat); % differences between true value and approximation of the integral
% exactsucc = mean(trueerrormat<=abstol,1); %percentage of successful instants
succnowarn = mean((trueerrormat<=abstol)&(~conechangemat),1); %percentage of successful instants for which the functions are in the cone
succwarn = mean((trueerrormat<=abstol)&(conechangemat),1);    %percentage of successful instants for which the functions are not in the cone
failnowarn = mean((trueerrormat>abstol)&(~conechangemat),1);  %percentage of failed instants for which the functions are in the cone
failwarn = mean((trueerrormat>abstol)&(conechangemat),1); %percentage of failed instants for which the functions are not in the cone
timeave = mean(timemat);
nave = mean(npointsmat);


%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
[fileID, fullPath] = gail.open_txt('WorkoutIntegralOutput', 'WorkoutIntegralTest');
fprintf(fileID,'\n          Success    Success   Failure    Failure    number      time');
fprintf(fileID,'\n ninit   No Warning  Warning  No Warning  Warning   of points   elasped\n');
for i=1:nninit
    fprintf(fileID,['%5.0f %10.2f%%' ...
        '%10.2f%% %8.2f%% %8.2f%% %9.2f %8.2f \n'],...
        [ninitvec(i) 100*[succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)] nave(i) timeave(i)]);
end
fclose(fileID);

%% Save Output

time = toc(tstart);

gail.save_mat('WorkoutIntegralOutput','WorkoutIntegralTest',true, nrep,time,...
        succnowarn,succwarn,failnowarn,failwarn,nave,timeave);
warning on
warning('on','GAIL:integral_g:exceedbudget');
warning('on','GAIL:integral_g:peaky');



end