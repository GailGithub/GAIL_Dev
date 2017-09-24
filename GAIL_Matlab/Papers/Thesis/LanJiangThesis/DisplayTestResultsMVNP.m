%This function is to plot the figures in Lan Jiang' thesis
% coloroption --- 'black' or 'color'
% Call the function as below:
% DisplayTestResultsMVNP('black') would plot
% example 1 using color marker.
function DisplayTestResultsMVNP(coloroption)
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
%load TestcubMCgFixedCovon-MVNP-N500d4abstol0.0001rel0.001-2015-09-17-23-38-21.mat
load TestcubMCgFixedCovon-MVNP-N500d14abstol0.0001rel0.001-2015-11-29-22-15-18.mat
%load TestcubMCon-MVNP-N500d6abstol0.0001rel0.0001-2015-09-08-14-30-41.mat

%load TestcubMCgFixedCovon-product-N500d2abstol0.001rel0-2015-11-24-16-01-55.mat

%load TestcubMCgFixedCovon-product-N500d17abstol0.001rel0.001-2015-11-26-00-53-24.mat
param.nrep = 500;
plotTest.plotcolor=coloroption;
plotTest.logerrtollo=-5;
plotTest.logerrtolhi=2;
plotTest.logtimelo=-4;
plotTest.logtimehi=3;
plotTest.errtollowlimit=10^plotTest.logerrtollo;
plotTest.errtolhilimit=10^plotTest.logerrtolhi;
plotTest.timelowlimit=10^plotTest.logtimelo;
plotTest.timehilimit=10^plotTest.logtimehi;
plotTest.linewidth=2;
plotTest.nrep=param.nrep;
plotTest.namepref=param.funtype;

%% Plot iid results
test.whichsample = 'iid';
if any(strcmp('iid',test.whichsample))
    plotTest.err=res.iiderr;
    plotTest.errtol = res.iiderr./res.tol;
    plotTest.time=res.iidtime;
    plotTest.exit=res.iidexit;
    %plotTest.kurtmax=res.kurtmax;
    %plotTest.kurtvec = res.kurtvec;
    plotTest.name=[plotTest.namepref 'iidErrTime'];
    plotTest.defaultcolor=[1,0,0];
    if any(strcmp('black',plotTest.plotcolor))
        plotTest.ptsize=150;
        plotTestcubMCgRELTOLblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
        plotTest.ptsize=150;
        plotTestcubMCgRELTOLcolor(plotTest,param)
    end
    %plotTest=rmfield(plotTest,'kmax');
end
%

%% Plot Sobol results
test.whichsample = 'cubsobol';
if any(strcmp('cubsobol',test.whichsample))
    plotTest.err=res.Sobolerr;
    plotTest.errtol = res.Sobolerr./res.tol;
    plotTest.time=res.Soboltime;
    plotTest.exit=res.Sobolexit;
    plotTest.name=[plotTest.namepref 'cubSobolErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
        plotTest.ptsize=150;
        plotTestcubMCgRELTOLblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
         plotTest.ptsize=200;
        plotTestcubMCgRELTOLcolor(plotTest,param)
    end
end

%% Plot Sobol heavy dutyresults
test.whichsample = 'cublattice';
if any(strcmp('cublattice',test.whichsample))
    plotTest.err=res.Latticeerr;
    plotTest.errtol = res.Latticeerr./res.tol;
    plotTest.time=res.Latticetime;
    plotTest.exit=res.Latticeexit;
    plotTest.name=[plotTest.namepref 'cubLatticeErrTime'];
    plotTest.defaultcolor=[1 0 0];
    if any(strcmp('black',plotTest.plotcolor))
        plotTest.ptsize=150;
        plotTestcubMCgRELTOLblack(plotTest,param)
    end
    if any(strcmp('color',plotTest.plotcolor))
         plotTest.ptsize=200;
        plotTestcubMCgRELTOLcolor(plotTest,param)
    end
end



end


