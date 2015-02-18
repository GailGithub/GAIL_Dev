%Run TestcubMC on the step function
clear all, close all
format compact

fun.funtype='geomean';
fun.S0=100;
fun.K=100;
fun.T=1;
fun.r=0.03;
param.measure='normal';
param.impyes=false;
param.abstol=2e-2;
param.reltol=0;
param.n_sigma=1e4;

test.nrep=500;
test.howoftenrep=10;
dimchoice=[1 2 4 8 16 32 64]';
ndim=size(dimchoice,1);
test.randch.dim=dimchoice(randi(ndim,test.nrep,1));
sigmin=0.1;
sigmax=0.7;
%test.randch.sigoverall=sigmin*(sigmax/sigmin).^rand(test.nrep,1);
test.randch.sigoverall=sigmin+(sigmax-sigmin).*rand(test.nrep,1);
test.randchoicefun=@randchoiceGeoCall;
test.whichsample={'cubLattice'};



%% TestCubatureDiffSettings
tempinitial=zeros(test.nrep,1);
res.dim=tempinitial;
if any(strcmp('cubLattice',test.whichsample))
    res.Latticeexit=tempinitial;
    res.LatticeQ=tempinitial;
    res.Latticeerr=tempinitial;
    res.Latticetime=tempinitial;
    res.Latticeneval=tempinitial;
end

for irep=1:test.nrep
    if round(irep/test.howoftenrep)==irep/test.howoftenrep, irep, end
    %keyboard
    param.sample='qmc';
    [testfunqmc,fun,param]= ...
          test.randchoicefun(fun,param,test.randch,irep);
    res.dim(irep)=param.dim;
    
    % Evaluate integral using cubLattice
        [q,out_param]=...
           cubLattice_g(testfunqmc,param.dim,...
           'abstol',param.abstol,'reltol',param.reltol,'measure',param.measure,'transform','Baker');
        %res.Latticeexit(irep)=out_param.overbudget;
        res.LatticeQ(irep)=q;
        res.Latticeerr(irep)=abs(param.exactintegral-q);
        res.Latticetime(irep)=out_param.time;
        res.Latticeneval(irep)=out_param.n;
end

timestamp=datestr(now,'yyyy-mm-dd-HH-MM');
save(['TestCubature-' fun.funtype '-' param.measure '-' timestamp ...
   '-N-' int2str(test.nrep) '-d-' int2str(param.dim)  ...
    '-tol-' num2str(param.abstol) '.mat'])


%% Display test results

%Display TestMCDiffSettings results

set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','tex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markersset(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
plotTest.logerrlo=-6;
plotTest.logerrhi=-1;

param.tol=param.abstol;
plotTest.plotcolor='color';
plotTest.logtimelo=-3;
plotTest.logtimehi=2;
plotTest.errlowlimit=10^plotTest.logerrlo;
plotTest.errhilimit=10^plotTest.logerrhi;
plotTest.timelowlimit=10^plotTest.logtimelo;
plotTest.timehilimit=10^plotTest.logtimehi;
plotTest.linewidth=2;
plotTest.nrep=test.nrep;
plotTest.namepref=fun.funtype;

% Plot Lattice results
if any(strcmp('cubLattice',test.whichsample))
    plotTest.err=res.Latticeerr;
    plotTest.time=res.Latticetime;
    plotTest.exit=res.Latticeexit;
    plotTest.name=[plotTest.namepref 'cubLatticeErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTest.ptsize=150;
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=400;
    plotTestcubMCcolor(plotTest,param)
    end
    Latticepercentright=mean(res.Latticeerr<=param.abstol)
end