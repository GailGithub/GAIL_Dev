%This script describes the Experiment 4:  
%  Functions with two local minimum points.
% 
%  Generates Table 3.4 in the thesis
%
%  Xin Tong, A Guaranteed, Adaptive, Automatic Algorithm for Univatiate
%  Function Minimization, July 2014.

%% Garbage collection and initialization
clearvars -except testCase  %clear all variables except testCase
close all 
tstart = tic;

%% Program parameters
TolXvec = [10^(-2) 10^(-4) 10^(-7)];
in_param.abstol = 0; %error tolerance
in_param.nmax = 10^7; %cost budget

%% Simulation parameters
nrep = 10000;
if (nrep >= 1000)
    warning('off','MATLAB:funmin_g:exceedbudget');
    warning('off','MATLAB:funmin_g:peaky');
end;
a1=5; b1=10; c1=0.5-0.5*rand(nrep,1);
a2=1; b2=10; c2=0.5+0.5*rand(nrep,1);
nTolX = length(TolXvec);

intnum = zeros(nrep,nTolX);
succfunmin = zeros(nrep,nTolX);
succfminbnd = zeros(nrep,nTolX);
xmin = zeros(nrep,nTolX);
exceedmat = zeros(nrep,nTolX);
exactsolu = zeros(nrep,nTolX); 

for i=1:nTolX
    in_param.TolX = TolXvec(i);
    for j=1:nrep
        f=@(x) -a1*exp(-(b1*(x-c1(j))).^2)-a2*exp(-(b2*(x-c2(j))).^2);
        exactsolu(j) = fminbnd(f,0,(c1(j)+c2(j))/2,optimset('TolX',1e-9));
        [fmin,out_param] = funmin_g(f,in_param);
        xmin(j,i)=fminbnd(f,0,1,optimset('TolX',in_param.TolX));
        intnum(j,i) = size(out_param.intervals,2);
        exceedmat(j,i) = out_param.exitflag;
        for k=1:intnum(j,i)
            if exactsolu(j) <= out_param.intervals(2,k) && exactsolu(j)... 
                    >= out_param.intervals(1,k) 
                succfunmin(j,i) = 1;
            end
        end
        if abs(xmin(j,i)-exactsolu(j)) <= in_param.TolX
            succfminbnd(j,i) = 1;
        end
    end
end

probfunmin=mean(succfunmin,1); %probability find the solution by funmin_g 
probnowarn=mean(succfunmin&(~exceedmat),1); 
probwarn=mean(succfunmin&(exceedmat),1); 
probfminbnd=mean(succfminbnd,1); 


%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('           Success    Success    Success    Success')
display(' TolX               No Warning   Warning    fminbnd' )
for i=1:nTolX
    display(sprintf([ '%1.0e    %7.2f%%   %7.2f%%   %7.2f%%   %7.2f%% '],...
        [TolXvec(i) 100*[probfunmin(i) probnowarn(i) probwarn(i)...
        probfminbnd(i)]])) 
end

%% Save Output
% [GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
% path = strcat(GAILPATH,'OutputFiles' , PATHNAMESEPARATOR);
% filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
%                   'TwoExtremeTest-',...
%                   datestr(now,'yyyymmdd'),'.mat');
% save(filename)

gail.save_mat('UniFunMinOutput', 'TwoExtremeTest',TolXvec,probfunmin,...
    probnowarn,probwarn,probfminbnd);

toc(tstart)

%% The following output was obtained on 2014-May
%            Success    Success    Success    Success
%  TolX               No Warning   Warning    fminbnd
% 1e-02     100.00%    100.00%      0.00%     67.28% 
% 1e-04     100.00%    100.00%      0.00%     67.28% 
% 1e-07     100.00%      0.00%    100.00%     67.28% 
% Elapsed time is 10873.821636 seconds.

