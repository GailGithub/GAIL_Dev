%% Keister's example of multidimensional integration
% 
% B. D. Keister, Multidimensional quadrature algorithms, _Computers in
% Physics_, *10*, pp. 119-122, 1996, presents the following
% multidimensional integral, inspired by a physics application:
%
% \[ I = \int_{\mathbb{R}^d} \cos(\lVert \boldsymbol{x} \rVert)
% \exp(-\lVert \boldsymbol{x} \rVert^2) \, \mathrm{d} \boldsymbol{x},
% \qquad d = 1, 2, \ldots. \]

%% Expressing the integral as an expectation
% Let's evaluate the integral using Monte Carlo cubature.  We first note
% that the change of variable \(\boldsymbol{t} = \boldsymbol{x}/a\)
% transforms this integral into
%
% \begin{align*} I &= \int_{\mathbb{R}^d} \cos(a \lVert \boldsymbol{t}
% \rVert) \exp(-a^2 \lVert \boldsymbol{t} \rVert^2) \, a^d \mathrm{d}
% \boldsymbol{t}, \qquad a > 0, \\ & = \int_{\mathbb{R}^d}
% \underbrace{(2\pi a^2)^{d/2} \cos(a \lVert \boldsymbol{t} \rVert)
% \exp((1/2-a^2) \lVert \boldsymbol{t} \rVert^2)}_{f(\boldsymbol{t})}
% \times \underbrace{\frac{\exp(-\lVert \boldsymbol{t} \rVert^2/2)}
% {(2\pi)^{d/2}}}_{\varrho(\boldsymbol{t})} \, \mathrm{d} \boldsymbol{t} \\
% & = \mathbb{E}[f(\boldsymbol{T})], \qquad \text{where } \boldsymbol{T} \sim \mathcal{N}(\boldsymbol{0},
% \mathsf{I}). \end{align*}

%% Evaluating the integral using |meanMC_g|
% To find \(I\) by Monte Carlo methods we define an anonymous function
% \(f\) as follows:

function [succTable,avgAbsErrTable,timeTable,nSampleTable,warnTable] = KeisterCubatureExampleWiley(nRep) 
%make it a function to not overwrite other variables

if nargin < 1
   nRep = 10;
end
gail.InitializeDisplay %initialize the display parameters
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
dvec = [3 8]; %vector of dimensions
abstol = [0.005 0.05]; %absolute error tolerance
reltol = [0 0]; %relative error tolerance
nd = length(dvec);
a = 1/sqrt(2); %"best" value of a

%%
% Next we call |meanMC_g| or |cubMC_g|
IMCvec(nRep,nd) = 0; %vector of answers
timeMC(nRep,nd) = 0;
nSampleMC(nRep,nd) = 0;
nWarnMC(nRep,nd) = 0;
tic
for i=1:nd
   d = dvec(i);
   for k = 1:nRep
      gail.print_iterations(k, 'k', true);
%      [IMCvec(k,i),out] = meanMC_g(@(n) f(randn(n,d),a,d),abstol(i),reltol(i));
      [IMCvec(k,i),out] = cubMC_g(@(x) f(x,a,d),[-inf(1,d); inf(1,d)], ...
         'normal',abstol(i),reltol(i));
       if out.exitflag > 0
          nWarnMC(k,i) = 1;
      end
      timeMC(k,i) = out.time;
      nSampleMC(k,i) = out.ntot;
   end
end
disp(' ')
toc
nSampleMC = mean(nSampleMC)
timeMC = mean(timeMC)

%% Now we do Sobol
ISobvec(nRep,nd) = 0; %vector of answers
timeSob(nRep,nd) = 0;
nSampleSob(nRep,nd) = 0;
nWarnSob(nRep,nd) = 0;
tic
for i=1:nd
   d = dvec(i);
   for k = 1:nRep
      gail.print_iterations(k, 'k', true);
      [ISobvec(k,i),out] = cubSobol_g(@(x) f(x,a,d),[-inf(1,d); inf(1,d)], ...
          'normal',abstol(i),reltol(i));
      if out.exitflag > 0
          nWarnSob(k,i) = 1;
      end
      timeSob(k,i) = out.time;
      nSampleSob(k,i) = out.n;
   end
end
disp(' ')
toc
timeSob = mean(timeSob)
nSampleSob = mean(nSampleSob)

%% Finally we we do lattice
ILatvec(nRep,nd) = 0; %vector of answers
timeLat(nRep,nd) = 0;
nSampleLat(nRep,nd) = 0;
nWarnLat(nRep,nd) = 0;
tic
for i=1:nd
    d = dvec(i);
    for k = 1:nRep
        gail.print_iterations(k, 'k', true);
        [ILatvec(k,i),out] = cubLattice_g(@(x) f(x,a,d),[-inf(1,d); inf(1,d)], ...
            'normal',abstol(i),reltol(i));
        if out.exitflag > 0
            nWarnLat(k,i) = 1;
        end
        timeLat(k,i) = out.time;
        nSampleLat(k,i) = out.n;
   end
end
disp(' ')
toc
timeLat = mean(timeLat)
nSampleLat = mean(nSampleLat)

%% Checking the real error
% There is a way to get the value of this integral to machine precision
% using the function |Keistertrue|
%
% <include>Keistertrue.m</include>

Ivec = repmat(Keistertrue(dvec),nRep,1);
absErrMC = abs(Ivec-IMCvec);
succMC = mean(absErrMC <= repmat(abstol,nRep,1))
avgAbsErrMC = mean(absErrMC)
absErrSob = abs(Ivec-ISobvec);
succSob = mean(absErrSob <= repmat(abstol,nRep,1))
avgAbsErrSob = mean(absErrSob)
absErrLat = abs(Ivec-ILatvec);
succLat = mean(absErrLat <= repmat(abstol,nRep,1))
avgAbsErrLat = mean(absErrLat)

%% Print number of warnings from each algorithms
warnMC = sum(nWarnMC,1);
warnLat = sum(nWarnLat,1);
warnSob = sum(nWarnSob,1);
disp(['warnings issued by cubMC_g: ', num2str(sum(nWarnMC,1))])
disp(['warnings issued by cubLattice_g: ', num2str(sum(nWarnLat,1))])
disp(['warnings issued by cubSobol_g: ', num2str(sum(nWarnSob,1))])


%outFileName = ['KeisterCubExWileyDataNRep' int2str(nRep) ...
%   datestr(now,'-yyyy-mm-dd-HH-MM-SS') '.mat'];
%save(outFileName)

outFileName = gail.save_mat('MC_StoppingCriteriaOutput',['KeisterCubExWileyDataNRep' int2str(nRep)],...
    true, abstol,...
    avgAbsErrMC, avgAbsErrLat, avgAbsErrSob, succMC, succLat, succSob,...
    timeMC, timeLat, timeSob, nSampleMC, nSampleLat, nSampleSob, ...
    warnMC, warnLat, warnSob);

KeisterCubExWileyOut(outFileName)

%% Prepare outputs
rowNames = {'cubMC_g','cubLattice_g','cubSobol_g'};
colNames = {'d_eq_3','d_eq_8'};
succTable = array2table([succMC; succLat; succSob],'RowNames',rowNames,'VariableNames',colNames);
%succTable(1,1)
avgAbsErrTable = array2table([avgAbsErrMC; avgAbsErrSob; avgAbsErrLat],'RowNames',rowNames,'VariableNames',colNames);
timeTable = array2table([timeMC; timeLat; timeSob],'RowNames',rowNames,'VariableNames',colNames);
nSampleTable = array2table([nSampleMC; nSampleLat; nSampleSob],'RowNames',rowNames,'VariableNames',colNames);
warnTable = array2table([warnMC; warnLat; warnSob],'RowNames',rowNames,'VariableNames',colNames);

%%
%
% _Author: Fred J. Hickernell_

   