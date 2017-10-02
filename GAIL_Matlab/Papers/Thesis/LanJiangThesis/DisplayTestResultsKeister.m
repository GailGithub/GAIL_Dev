%This function is to plot the figures in Lan Jiang' thesis
% coloroption --- 'black' or 'color'
% Call the function as below:
% DisplayTestResultsKeister('black') would plot
% example 1 using color marker.
function DisplayTestResultsKeister(coloroption)
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
load TestcubMCon-Keister-uniform-N500d18abstol0.001rel0.001-2015-08-24-21-47-53.mat
plotTest.plotcolor=coloroption;
plotTest.logerrtollo=-5;
plotTest.logerrtolhi=4;
plotTest.logtimelo=-4;
plotTest.logtimehi=3;
plotTest.errtollowlimit=10^plotTest.logerrtollo;
plotTest.errtolhilimit=10^plotTest.logerrtolhi;
plotTest.timelowlimit=10^plotTest.logtimelo;
plotTest.timehilimit=10^plotTest.logtimehi;
plotTest.linewidth=2;
plotTest.nrep=test.nrep;
plotTest.namepref='gaussian';
%% Plot iid results
test.whichsample = 'iid';
if any(strcmp('iid',test.whichsample))
    plotTest.err=res.iiderr;
    plotTest.errtol = res.iiderr./res.tol;
    plotTest.time=res.iidtime;
    plotTest.exit=res.iidexit;
    plotTest.kurtmax=res.iidkurtmax;
    %plotTest.kurtvec = res.exactkurtosis;
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
    plotTest=rmfield(plotTest,'kurtmax');
end
%

%% Plot cubSobol_G results
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

%% Plot cubLattice_g dutyresults
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

% %% Plot cubSobol_G results
% test.whichsample = 'integral';
% if any(strcmp('integral',test.whichsample))
%     plotTest.err=res.integralerr;
%     plotTest.errtol = res.integralerr./res.tol;
%     plotTest.time=res.integraltime;
%     plotTest.name=[plotTest.namepref 'integralErrTime'];
%     plotTest.defaultcolor=[1 0 0];
%     if any(strcmp('black',plotTest.plotcolor))
%         plotTest.ptsize=150;
%         plotTestcubMCgRELTOLblack(plotTest,param)
%     end
%     if any(strcmp('color',plotTest.plotcolor))
%          plotTest.ptsize=200;
%         plotTestcubMCgRELTOLcolor(plotTest,param)
%     end
% end
% 
% %% Plot chebfun results
% test.whichsample = 'chebfun';
% if any(strcmp('chebfun',test.whichsample))
%     plotTest.err=res.chebfunerr;
%     plotTest.errtol = res.chebfunerr./res.tol;
%     plotTest.time=res.chebfuntime;
%     plotTest.name=[plotTest.namepref 'chebfunErrTime'];
%     plotTest.defaultcolor=[1 0 0];
%     if any(strcmp('black',plotTest.plotcolor))
%         plotTest.ptsize=150;
%         plotTestcubMCgRELTOLblack(plotTest,param)
%     end
%     if any(strcmp('color',plotTest.plotcolor))
%          plotTest.ptsize=200;
%         plotTestcubMCgRELTOLcolor(plotTest,param)
%     end
% end


end


