%Display TestMCDiffSettings results
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markers
%load TestCubature-Genz-uniform-2014-08-06-17-15-N-1000-d-1-tol-0.001.mat
%load TestCubature-Genz-uniform-2014-08-07-08-43-N-1000-d-7-tol-0.001.mat %mod
%load TestCubature-Genz-uniform-2014-08-07-08-45-N-1000-d-6-tol-0.001.mat %mod
%load TestCubature-Genz-uniform-2014-08-07-09-40-N-1000-d-5-tol-0.001.mat
%load TestCubature-Genz-uniform-2014-08-07-09-42-N-1000-d-5-tol-0.001.mat
%load TestCubature-Genz-uniform-2014-08-07-16-34-N-1000-d-7-tol-0.001.mat
%load TestCubature-Genz-uniform-2014-08-07-16-40-N-1000-d-7-tol-0.001.mat
%load TestCubature-Genz-uniform-2014-08-07-16-48-N-1000-d-1-tol-0.001.mat
%load TestCubature-Genz-uniform-2014-08-07-17-08-N-1000-d-10-tol-0.001.mat
%load TestCubature-Keister-uniform-2014-08-08-16-30-N-1000-d-7-tol-0.001.mat
%load TestCubature-Keister-uniform-2014-08-08-17-03-N-1000-d-19-10a-13-tol-0.001.mat
%load TestCubature-Keister-uniform-2014-08-08-18-05-N-1000-d-19-10a-10-tol-0.001.mat
%load TestCubature-Keister-uniform-2014-08-08-18-08-N-1000-d-19-10a-14-tol-0.001.mat
load TestCubature-Keister-uniform-2014-08-12-10-50-N-1000-d-19-10a-14-tol-0.001.mat
plotTest.plotcolor='black';
plotTest.logerrlo=-6;
plotTest.logerrhi=-1;
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

%% Plot iid results
if any(strcmp('iid',test.whichsample))
    plotTest.err=res.iiderr;
    plotTest.time=res.iidtime;
    plotTest.exit=res.iidexit;
    plotTest.kurtmax=res.iidkurtmax;
    plotTest.name=[plotTest.namepref 'iidErrTime'];
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
end
% 
%% Plot iid heavy duty results
if any(strcmp('iidheavy',test.whichsample))
    plotTest.err=res.iidheavyerr;
    plotTest.time=res.iidheavytime;
    plotTest.exit=res.iidheavyexit;
    plotTest.kurtmax=res.iidheavykurtmax;
    plotTest.name=[plotTest.namepref 'iidheavyErrTime'];
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
end

%% Plot Sobol results
if any(strcmp('cubSobol',test.whichsample))
    plotTest.err=res.Sobolerr;
    plotTest.time=res.Soboltime;
    plotTest.exit=res.Sobolexit;
    plotTest.name=[plotTest.namepref 'cubSobolErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTest.ptsize=200;
    plotTestcubMCblack_FJH(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTest.ptsize=400;
    plotTestcubMCcolor(plotTest,param)
    end
    Sobolsuccess=mean(res.Sobolerr<=param.tol)
end

