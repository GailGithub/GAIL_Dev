%Display TestMCDiffSettings results
clearvars %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',16,'defaulttextfontsize',16) %make font larger
set(0,'defaultLineLineWidth',4) %thick lines
set(0,'defaultTextInterpreter','tex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markersset(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
% load TestCubature-geomean-normal-2015-06-26-13-15-N-500-d-4-tol-0.02.mat
% plotTest.logerrlo=-5;
% plotTest.logerrhi=0;

load TestCubature-Keister-uniform-N-500-cubLattice_rel.mat
% load TestCubature-Keister-uniform-N-500-cubSobol_rel.mat
plotTest.print = 0;
plotTest.logerrhi=3;
plotTest.logerrlo=-4;
param.tol=param.abstol;
plotTest.plotcolor='color';
plotTest.logtimelo=-3;
plotTest.logtimehi=0;
plotTest.errlowlimit=10^plotTest.logerrlo;
plotTest.errhilimit=10^plotTest.logerrhi;
plotTest.timelowlimit=10^plotTest.logtimelo;
plotTest.timehilimit=10^plotTest.logtimehi;
plotTest.linewidth=2;
plotTest.nrep=test.nrep;
plotTest.namepref=fun.funtype;
if strcmp(fun.funtype,'step') 
plotTest.kurtvec=res.exactkurtosis;
plotTest.namepref=[plotTest.namepref 'd=' int2str(param.dim)];
end
if strcmp(fun.funtype,'gaussian') 
plotTest.kurtvec=res.exactkurtosis;
plotTest.namepref=[plotTest.namepref 'd=' int2str(param.dim)];
end

%% Plot Sobol results
if any(strcmp('cubSobol',test.whichsample))
    plotTest.err=res.Sobolerr./gail.tolfun(param.abstol,param.reltol,param.theta,res.Sobolexact,param.toltype);;
    plotTest.time=res.Soboltime;
    plotTest.exit=res.Sobolexit;
    plotTest.name=[plotTest.namepref 'cubSobolErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
        error('Not plot in BW')
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=200;
    plotTestColor_rel(plotTest,param)
    end
    Sobolpercentright=mean(plotTest.err<=1)
end

%% Plot Lattice results
if any(strcmp('cubLattice',test.whichsample))
    plotTest.err=res.Latticeerr./gail.tolfun(param.abstol,param.reltol,param.theta,res.Latticeexact,param.toltype);
    plotTest.time=res.Latticetime;
    plotTest.exit=res.Latticeexit;
    plotTest.name=[plotTest.namepref 'cubLatticeErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
        error('Not plot in BW')
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=200;
    plotTestColor_rel(plotTest,param)
    end
    Latticepercentright=mean(plotTest.err<=1)
end

