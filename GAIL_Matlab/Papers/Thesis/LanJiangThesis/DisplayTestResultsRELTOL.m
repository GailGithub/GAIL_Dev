%This function is to plot the figures in Lan Jiang' thesis
% coloroption --- 'black' or 'color'
% Call the function as below:
% DisplayTestResultsRELTOL({'ex1'},'black') would plot
% example 1 using color marker.

function DisplayTestResultsRELTOL(NumericalExample,coloroption)
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
for i = 1:length(NumericalExample)
    switch NumericalExample{i}
        case 'ex1'
            load TestcubMCon-gaussian-uniform-N500d1abstol0.001-2015-11-27-18-31-04.mat
            %load TestcubMCon-gaussian-uniform-N500d1abstol0.001-2015-11-16-12-34-44.mat
            %load TestcubMCon-MVNP-N500d4abstol0.0001rel0.0001-2015-09-11-22-10-44.mat
            %load TestcubMCon-MVNP-N5d5abstol0.0001rel0.0001-2015-09-10-01-20-12.mat
            %load TestcubMCon-Keister-uniform-N500d18abstol0.001rel0.001-2015-08-24-21-47-53.mat
            %load TestcubMCon-Keister-uniform-N500d17abstol0.001rel0.01-2015-08-22-02-18-46.mat
            %load TestcubMCon-gaussianker-uniform-N500d3abstol1e-08rel0.001-2015-08-15-23-40-57.mat
            %load TestcubMCon-gaussianker-uniform-N500d1abstol1e-08rel0.0001-2015-08-18-23-31-58.mat
            %load TestcubMCon-gaussianker-uniform-N500d1abstol1e-08rel0.001-2015-08-19-00-12-52.mat
        case 'ex2'
            if exist('TestcubMCon-gaussian-uniform-Out-17-Aug-2012_17.46.40N500d6tol0.001.mat')
                load TestcubMCon-gaussian-uniform-Out-17-Aug-2012_17.46.40N500d6tol0.001.mat
                % this MAT file is to plot figure 5 in paper
            else
                warning(['TestcubMCon-gaussian-uniform-Out-17-Aug-2012_17.46.40N500d6tol0.001.mat does not exist. '...
                    'need to call function RunTestcubMConGaussian to produce the MAT file.'])
                [res,test,fun,param] = RunTestcubMConGaussian;
            end
        case 'ex3'
            if exist('TestcubMCon-geomean-Out-17-Aug-2012_20.38.24N500d1tol0.05.mat')
                load TestcubMCon-geomean-Out-17-Aug-2012_20.38.24N500d1tol0.05.mat
           % this MAT file is to plot figure 6 in paper
            else
                warning(['TestcubMCon-geomean-Out-17-Aug-2012_20.38.24N500d1tol0.05.mat does not exist. '...
                    'need to call function RunTestcubMConGeoAsianCall to produce the MAT file.'])
                [res,test,fun,param] = RunTestcubMConGeoAsianCall;
            end
            
    end    
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
    plotTest.nrep=test.nrep;
    plotTest.namepref=fun.funtype;
    if strcmp(fun.funtype,'gaussianker')
        plotTest.kurtvec=res.exactkurtosis;
        plotTest.namepref=[plotTest.namepref 'd=' int2str(param.dim)];
    end
    
    %% Plot iid results
    if any(strcmp('iid',test.whichsample))
        plotTest.err=res.iiderr;
        plotTest.errtol = res.iiderr./res.tol;
        plotTest.time=res.iidtime;
        plotTest.exit=res.iidexit;
        plotTest.kurtmax=res.iidkurtmax;
        plotTest.name=[plotTest.namepref 'iidErrTime'];
        plotTest.defaultcolor=[1,0,0];
        if any(strcmp('black',plotTest.plotcolor))
            plotTest.ptsize=150;
            plotTestcubMCgRELTOLblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTest.ptsize=400;
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
        plotTest=rmfield(plotTest,'kurtmax');
    end
    %
    %% Plot iid heavy duty results
    if any(strcmp('iidheavy',test.whichsample))
        plotTest.err=res.iidheavyerr;
        plotTest.errtol = res.iidheavyerr./res.tol;
        plotTest.time=res.iidheavytime;
        plotTest.exit=res.iidheavyexit;
        plotTest.kurtmax=res.iidheavykurtmax;
        plotTest.name=[plotTest.namepref 'iidheavyErrTime'];
        if any(strcmp('black',plotTest.plotcolor))
            plotTestcubMCgRELTOLblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
        plotTest=rmfield(plotTest,'kurtmax');
    end
    
    %% Plot Sobol results
    if any(strcmp('cubSobol',test.whichsample))
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
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
    end
    
    %% Plot Sobol heavy dutyresults
    if any(strcmp('cubLattice',test.whichsample))
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
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
    end
    
    
    %% Plot chebfun and chebfunHeavy results
    if any(strcmp('chebfun',test.whichsample))
        plotTest.err=res.chebfunerr;
        plotTest.time=res.chebfuntime;
        plotTest.name=[plotTest.namepref 'chebfunErrTime'];
        plotTest.defaultcolor=[1 0 0];
        if any(strcmp('black',plotTest.plotcolor))
            plotTestcubMCgRELTOLblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
    end
    
    if any(strcmp('chebfunheavy',test.whichsample))
        plotTest.err=res.chebfunheavyerr;
        plotTest.time=res.chebfunheavytime;
        plotTest.name=[plotTest.namepref 'chebfunheavyErrTime'];
        plotTest.defaultcolor=[1 0 0];
        if any(strcmp('black',plotTest.plotcolor))
            plotTestcubMCgRELTOLblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
    end
        %% Plot chebfun and chebfunHeavy results
    if any(strcmp('integral',test.whichsample))
        plotTest.err=res.integralerr;
        plotTest.errtol = res.integralerr./res.tol;
        plotTest.time=res.integraltime;
        plotTest.name=[plotTest.namepref 'integralErrTime'];
        plotTest.defaultcolor=[1 0 0];
        if any(strcmp('black',plotTest.plotcolor))
            plotTestcubMCgRELTOLblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTestcubMCgRELTOLcolor(plotTest,param)
        end
    end
end
end


