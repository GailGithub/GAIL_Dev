%% Experiment 2: Bump test functions with abstol=0 & TolX=10^(-6)

%% Garbage collection and initialization
format long e 
clear all 
close all 
tstart = tic;

%% Program parameters
in_param.abstol = 0; %error tolerance
in_param.TolX = 10^(-6);
in_param.nmax = 10^7; %cost budget

%% Simulation parameters
nrep = 10000;
if (nrep >= 1000)
    warning('off','MATLAB:funmin_g:exceedbudget');
    warning('off','MATLAB:funmin_g:peaky');
end;
a = 10.^(-4+3*rand(nrep,1));
z = 2.*a+(1-4*a).*rand(nrep,1);
x0 = z-2*a;
x1 = z+2*a;
tauvec = [11 101 1001]; %cone condition tau
ntau = length(tauvec);
ratio = 1./a;
gnorm = 1./a;
exactmin = -1;
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
        [fmin,out_param] = funmin_g(f,in_param);
        ntrapmat(j,i) = out_param.npoints;
        newtaumat(j,i) = out_param.tau;
        estmin = fmin;
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

probinit = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); 
probfinl = mean(repmat(ratio,1,ntau)<=newtaumat,1); 
succnowarn = mean((truesolumat)&(~exceedmat),1); 
succwarn = mean((truesolumat)&(exceedmat),1);    
failnowarn = mean((~truesolumat)&(~exceedmat),1);   
failwarn = mean((~truesolumat)&(exceedmat),1);


%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:ntau
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%%'],...
        [tauvec(i) 100*[probinit(i) probfinl(i) succnowarn(i) ...
        succwarn(i) failnowarn(i) failwarn(i) ]])) 
end

%% Save Output
[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
path = strcat(GAILPATH,'OutputFiles' , PATHNAMESEPARATOR);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
                  'XToleranceTest-',...
                  datestr(now,'yyyymmdd'),'.mat');
save(filename)

toc(tstart)

%% The following output was obtained on 2014-May
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.52%->20.92%   20.92%      0.00%   79.08%    0.00%
%   101 32.80%->52.19%   52.19%      0.00%   47.81%    0.00%
%  1001 65.93%->84.92%   80.31%      4.64%   15.05%    0.00%
% Elapsed time is 6127.763619 seconds.


