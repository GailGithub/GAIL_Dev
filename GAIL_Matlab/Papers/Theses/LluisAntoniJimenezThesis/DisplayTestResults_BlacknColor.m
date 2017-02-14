%Display TestMCDiffSettings results
clearvars %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',16,'defaulttextfontsize',16) %make font larger
set(0,'defaultLineLineWidth',4) %thick lines
set(0,'defaultTextInterpreter','tex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markersset(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

% plotTest.logerrhi=-1;
% plotTest.logerrlo=-6;
% plotTest.logtimelo=-3;
% plotTest.logtimehi=0;
% load TestCubature-Keister-uniform-N-500-cubLattice.mat
% load TestCubature-Keister-uniform-N-500-cubLatticeCLT.mat
% load TestCubature-Keister-uniform-N-500-cubSobol.mat
% load TestCubature-Keister-uniform-N-500-cubSobolCLT.mat

plotTest.logerrlo=-5;
plotTest.logerrhi=0;
plotTest.logtimelo=-2;
plotTest.logtimehi=1;
% load TestCubature-geomean-normal-N-500-cubSobol.mat
% load TestCubature-geomean-normal-N-500-cubSobolCLT.mat
load TestCubature-geomean-normal-N-500-cubLattice.mat
% load TestCubature-geomean-normal-N-500-cubLatticeCLT.mat

plotTest.print = 0;
param.tol=param.abstol;
plotTest.plotcolor='color';

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
    plotTest.err=res.Sobolerr;
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
    plotTestColor(plotTest,param)
    end
    Sobolpercentright=mean(res.Sobolerr<=param.abstol)
end

%% Plot Lattice results
if any(strcmp('cubLattice',test.whichsample))
    plotTest.err=res.Latticeerr;
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
    plotTestColor(plotTest,param)
    end
    Latticepercentright=mean(res.Latticeerr<=param.abstol)
end

if any(strcmp('cubLatticeCLT',test.whichsample))
    plotTest.err=res.Latticeerr;
    plotTest.time=res.Latticetime;
    plotTest.exit=res.Latticeexit;
    plotTest.name=[plotTest.namepref 'cubLatticeCLTErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
        error('Not plot in BW')
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=200;
    plotTestColor(plotTest,param)
    end
    Latticepercentright=mean(res.Latticeerr<=param.abstol)
end

if any(strcmp('cubSobolCLT',test.whichsample))
    plotTest.err=res.Sobolerr;
    plotTest.time=res.Soboltime;
    plotTest.exit=res.Sobolexit;
    plotTest.name=[plotTest.namepref 'cubSobolCLTErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
        error('Not plot in BW')
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=200;
    plotTestColor(plotTest,param)
    end
    Sobolpercentright=mean(res.Sobolerr<=param.abstol)
end

