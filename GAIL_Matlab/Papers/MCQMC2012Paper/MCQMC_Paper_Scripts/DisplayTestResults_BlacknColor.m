%Display TestMCDiffSettings results
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
close all
clear all
load TestcubMCon-gaussian-uniform-Out-17-Aug-2012_13.12.36N500d1tol0.001.mat
%load TestcubMCon-gaussian-uniform-Out-01-Sep-2013_23.04.44N5d1tol0.001.mat
plotTest.plotcolor='color';
plotTest.logerrlo=-5;
plotTest.logerrhi=0;
plotTest.logtimelo=-3;
plotTest.logtimehi=3;
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
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)
    end
    plotTest=rmfield(plotTest,'kurtmax');
end

%% Plot Sobol results
if any(strcmp('Sobol',test.whichsample))
    plotTest.err=res.Sobolerr;
    plotTest.time=res.Soboltime;
    plotTest.exit=res.Sobolexit;
    plotTest.name=[plotTest.namepref 'SobolErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)
    end
end

%% Plot Sobol heavy dutyresults
if any(strcmp('Sobolheavy',test.whichsample))
    plotTest.err=res.Sobolheavyerr;
    plotTest.time=res.Sobolheavytime;
    plotTest.exit=res.Sobolheavyexit;
    plotTest.name=[plotTest.namepref 'SobolheavyErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)   
    end
end

%% Plot quad or quadgk results
if any(strcmp('quad',test.whichsample))
    plotTest.err=res.quaderr;
    plotTest.time=res.quadtime;
    plotTest.name=[plotTest.namepref 'quadErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)
    end
end

if any(strcmp('quadgk',test.whichsample))
    plotTest.err=res.quadgkerr;
    plotTest.time=res.quadgktime;
    plotTest.name=[plotTest.namepref 'quadgkErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)
    end
end

%% Plot chebfun and chebfunHeavy results
if any(strcmp('chebfun',test.whichsample))
    plotTest.err=res.chebfunerr;
    plotTest.time=res.chebfuntime;
    plotTest.name=[plotTest.namepref 'chebfunErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)
    end
end

if any(strcmp('chebfunheavy',test.whichsample))
    plotTest.err=res.chebfunheavyerr;
    plotTest.time=res.chebfunheavytime;
    plotTest.name=[plotTest.namepref 'chebfunheavyErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
    plotTestcubMCblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
    plotTestcubMCcolor(plotTest,param)
    end
end


