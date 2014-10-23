%% Experiment 1: Bump test functions with abstol=10^(-8) & TolX=0

%% Garbage collection and initialization
format long e 
clear all 
close all 
tstart = tic;

%% Program parameters
in_param.abstol = 10^(-8); %error tolerance
in_param.TolX = 0;
in_param.nmax = 10^7; %cost budget

%% Simulation parameters
nrep = 10000;
if (nrep >= 1000)
    warning('off','MATLAB:funmin_g:exceedbudget');
    warning('off','MATLAB:funmin_g:peaky');
end;
a = 10.^(-4+3*rand(nrep,1));
z = 2.*a+(1-4*a).*rand(nrep,1);
tauvec = [11 101 1001]; %cone condition tau
ntau = length(tauvec);
ratio = 1./a;
gnorm = 1./a;
exactmin = -1;

%% Simulation
ntrapmat = zeros(nrep,ntau);
trueerrormat = ntrapmat;
truesolumat = ntrapmat;
newtaumat = ntrapmat;
tauchangemat = ntrapmat;
exceedmat = ntrapmat;

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
        exceedmat(j,i) = out_param.exceedbudget;
    end
end

probinit = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); 
probfinl = mean(repmat(ratio,1,ntau)<=newtaumat,1); 
succnowarn=mean((trueerrormat<=in_param.abstol|truesolumat)&(~exceedmat),1); 
succwarn=mean((trueerrormat<=in_param.abstol|truesolumat)&(exceedmat),1);    
failnowarn=mean((trueerrormat>in_param.abstol&~truesolumat)&(~exceedmat),1);  
failwarn=mean((trueerrormat>in_param.abstol&~truesolumat)&(exceedmat),1);  

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:ntau
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%% '],...
        [tauvec(i) 100*[probinit(i) probfinl(i) succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i)]])) 
end

 
%% Save Output
[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
path = strcat(GAILPATH,'OutputFiles' , PATHNAMESEPARATOR);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
                  'ErrorToleranceTest-',...
                  datestr(now,'yyyymmdd'),'.mat');
save(filename)

toc(tstart)

%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.24%->57.53%   47.82%      3.87%   42.47%    5.84% 
%   101 33.03%->75.35%   66.50%      3.32%   24.65%    5.53% 
%  1001 66.70%->93.64%   82.23%      5.25%    6.36%    6.16% 

