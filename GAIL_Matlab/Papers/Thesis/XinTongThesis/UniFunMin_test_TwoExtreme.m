%This function describes the Experiment 4
% 
%  Generates Table 3.4 in the thesis with TolX=[10^(-2) 10^(-4) 10^(-7)], 
%  nrep=10000 and nmax=10^7 
%
%  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for Univariate
%  Function Minimization. MS thesis, Illinois Institute of Technology,
%  2014.

function [TolXvec,prob]=UniFunMin_test_TwoExtreme(nrep,TolX,nmax)

%% Program parameters
TolXvec = TolX;
in_param.abstol = 0; %error tolerance
in_param.nmax = nmax; %cost budget
tstart = tic;

%% Simulation parameters
n = nrep;
if (n >= 100)
    warning('off','GAIL:funmin01_g:exceedbudget');
    warning('off','GAIL:funmin01_g:peaky');
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
        [fmin,out_param] = funmin01_g(f,in_param);
        xmin(j,i)=fminbnd(f,0,1,optimset('TolX',in_param.TolX));
        intnum(j,i) = size(out_param.intervals,2);
        exceedmat(j,i) = out_param.exceedbudget;
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

prob.probfunmin=mean(succfunmin,1); %probability find the solution by funmin_g 
prob.probnowarn=mean(succfunmin&(~exceedmat),1); 
prob.probwarn=mean(succfunmin&(exceedmat),1); 
prob.probfminbnd=mean(succfminbnd,1); 

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('           Success    Success    Success    Success')
display(' TolX               No Warning   Warning    fminbnd' )
for i=1:nTolX
    display(sprintf([ '%1.0e    %7.2f%%   %7.2f%%   %7.2f%%   %7.2f%% '],...
        [TolXvec(i) 100*[prob.probfunmin(i) prob.probnowarn(i) ...
        prob.probwarn(i) prob.probfminbnd(i)]])) 
end

%% Save output
gail.save_mat('UniFunMinOutput', 'TwoExtremeTest',true,TolXvec,prob,nTolX);

toc(tstart)

warning('on','GAIL:funmin01_g:exceedbudget');
warning('on','GAIL:funmin01_g:peaky');

%% The following output was obtained on 2015-Feb
%            Success    Success    Success    Success
%  TolX               No Warning   Warning    fminbnd
% 1e-02     100.00%    100.00%      0.00%     66.76% 
% 1e-04     100.00%    100.00%      0.00%     66.76% 
% 1e-07     100.00%      0.00%    100.00%     66.76% 
