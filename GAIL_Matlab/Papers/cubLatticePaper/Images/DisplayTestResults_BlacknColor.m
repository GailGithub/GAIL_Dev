%Display TestMCDiffSettings results
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','tex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markersset(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
%load TestCubature-Keister-uniform-2014-09-12-15-46-N-1000-d-19-tol-0.001.mat
%load TestCubature-geomean-normal-2014-09-12-16-04-N-1000-d-1-tol-0.002.mat
%load TestCubature-Keister-uniform-2014-09-25-18-20-N-100-d-3-tol-0.001.mat
%load TestCubature-geomean-normal-2014-10-01-13-19-N-500-d-32-tol-0.02.mat %IID results
load TestCubature-geomean-normal-2014-10-16-11-25-N-500-d-4-tol-0.02.mat % Results for the paper
% plotTest.logerrlo=-5;
% plotTest.logerrhi=0;
%load TestCubature-Keister-uniform-2014-10-01-20-34-N-500-d-2-tol-0.002.mat %IID results
%load TestCubature-Keister-uniform-2014-10-02-16-17-N-500-d-9-tol-0.002.mat %Sobol and lattice results d<=19
%load TestCubature-Keister-uniform-2014-10-02-16-16-N-500-d-1-tol-0.002.mat %Sobol and lattice results d<=4
%load TestCubature-Keister-uniform-2014-10-06-14-17-N-500-d-1-tol-0.002.mat %Show Tony cubLattice
%load TestCubature-Keister-uniform-2014-10-06-14-24-N-500-d-12-tol-0.002.mat %Experiment with fudge
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
if strcmp(fun.funtype,'step') 
plotTest.kurtvec=res.exactkurtosis;
plotTest.namepref=[plotTest.namepref 'd=' int2str(param.dim)];
end
if strcmp(fun.funtype,'gaussian') 
plotTest.kurtvec=res.exactkurtosis;
plotTest.namepref=[plotTest.namepref 'd=' int2str(param.dim)];
end

%% Plot iid results
if any(strcmp('iid',test.whichsample))
    plotTest.err=res.iiderr;
    plotTest.time=res.iidtime;
    plotTest.exit=res.iidexit;
    plotTest.kurtmax=res.iidkurtmax;
    plotTest.name=[plotTest.namepref 'iidErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1,0,0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTest.ptsize=150;
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=400;
    plotTestcubMCcolor(plotTest,param)
    end
    plotTest=rmfield(plotTest,'kurtmax');
    IIDpercentright=mean(res.iiderr<=param.abstol)
end
% 
%% Plot iid heavy duty results
if any(strcmp('iidheavy',test.whichsample))
    plotTest.err=res.iidheavyerr;
    plotTest.time=res.iidheavytime;
    plotTest.exit=res.iidheavyexit;
    plotTest.kurtmax=res.iidheavykurtmax;
    plotTest.name=[plotTest.namepref 'iidheavyErrTime_d_' ...
       int2str(max(res.dim))];
    plotTest.defaultcolor=[1,0,0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTest.ptsize=150;
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=400;
    plotTestcubMCcolor(plotTest,param)
    end
    plotTest=rmfield(plotTest,'kurtmax');
    IIDHeavypercentright=mean(res.iidheavyerr<=param.tol)
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
    plotTest.ptsize=150;
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=400;
    plotTestcubMCcolor(plotTest,param)
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
    plotTest.ptsize=150;
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=400;
    plotTestcubMCcolor(plotTest,param)
    end
    Latticepercentright=mean(res.Latticeerr<=param.abstol)
end

