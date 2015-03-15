function [tauvec,prob] =workout_ErrXToleranceTest(nrep,abstol,TolX,nmax)
% user can choose number of iteration, absolut error tolerance, X
% tolerance, and cost budget nmax. 
%
% Experiment 3: Bump test functions with abstol=10^(-8), TolX=10^(-6), 
% nrep=10000 and nmax=10^7 


%% Program parameters
in_param.abstol = abstol; %error tolerance
in_param.TolX = TolX;
in_param.nmax = nmax; %cost budget

%% Simulation parameters
n = nrep;
if (n >= 100)
    warning('off','MATLAB:funmin_g:exceedbudget');
    warning('off','MATLAB:funmin_g:peaky');
end;
a = 10.^(-4+3*rand(n,1));
z = 2.*a+(1-4*a).*rand(n,1);
tauvec = [11 101 1001]; %cone condition tau
ntau = length(tauvec);
ratio = 1./a;
exactmin = -1;
exactsolu = z;

%% Simulation
ntrapmat = zeros(nrep,ntau);
trueerrormat = ntrapmat;
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
        in_param.nlo = (tauvec(i)+1)/2+1;
        in_param.nhi = in_param.nlo;
        [fmin,out_param] = funmin_g(f,in_param);
        ntrapmat(j,i) = out_param.npoints;
        newtaumat(j,i) = out_param.tau;
        estmin = fmin;
        trueerrormat(j,i) = abs(estmin-exactmin);
        tauchangemat(j,i) = out_param.tauchange;
        exceedmat(j,i) = out_param.exitflag;
        intnum(j,i) = size(out_param.intervals,2);
        for k=1:intnum(j,i)
            if exactsolu(j) <= out_param.intervals(2,k) && exactsolu(j) ... 
                    >= out_param.intervals(1,k) 
                truesolumat(j,i) = 1;
            end
        end
    end
end

warning('on','MATLAB:funmin_g:exceedbudget');
warning('on','MATLAB:funmin_g:peaky');

prob.probinit = mean(repmat(ratio,1,ntau)<=repmat(tauvec,nrep,1),1); 
prob.probfinl = mean(repmat(ratio,1,ntau)<=newtaumat,1); 
prob.succnowarn=mean((trueerrormat<=in_param.abstol|truesolumat)&(~exceedmat),1); 
prob.succwarn=mean((trueerrormat<=in_param.abstol|truesolumat)&(exceedmat),1);    
prob.failnowarn=mean((trueerrormat>in_param.abstol&~truesolumat)&(~exceedmat),1);  
prob.failwarn=mean((trueerrormat>in_param.abstol&~truesolumat)&(exceedmat),1);  

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('        Probability    Success   Success   Failure  Failure')
display(' tau      In Cone    No Warning  Warning No Warning Warning')
for i=1:ntau
    display(sprintf(['%5.0f %5.2f%%->%5.2f%% %7.2f%%' ...
        '%10.2f%% %7.2f%% %7.2f%% '],...
        [tauvec(i) 100*[prob.probinit(i) prob.probfinl(i) ...
        prob.succnowarn(i) prob.succwarn(i) prob.failnowarn(i)... 
        prob.failwarn(i)]])) 
end


%% Save output
gail.save_mat('WorkoutFunminOutput', 'workout_ErrXToleranceTest',tauvec,prob,ntau);

end
