%This function is to plot the figure 3-6 in MCQMC 2012 Paper
% Numerical Example --- three numerical examples {'ex1' 'ex2' 'ex3'} 
% coloroption --- 'black' or 'color' 
% Call the function as below:
% DisplayTestResults_BlacknColor({'ex1' 'ex2' 'ex3'},'black') would plot
% example 1-3 using black marker.
%
% DisplayTestResults_BlacknColor({'ex1'},'black') would plot
% example 1 using color marker.
%

function DisplayTestResults_BlacknColor(NumericalExample,coloroption)
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
for i = 1:length(NumericalExample)
    switch NumericalExample{i}
        case 'ex1'
            load TestcubMCon-gaussian-uniform-N500d1abstol0.001-2015-11-16-12-34-44.mat
            %load TestcubMCon-geomean-normal-N500d16abstol0.05-2015-09-13-15-36-31.mat
            %load TestcubMCon-Keister-uniform-N500d1abstol0.01rel0-2015-09-12-04-22-55.mat
            %load TestcubMCon-gaussian-uniform-N500d6abstol0.001-2015-08-11-21-13-18.mat
            %load TestcubMCon-gaussian-uniform-N500d1abstol0.001-2015-08-06-13-04-40.mat
%             if exist('TestcubMCon-gaussian-uniform-Out-17-Aug-2012_13.12.36N500d1tol0.001.mat')
%                 load TestcubMCon-gaussian-uniform-Out-17-Aug-2012_13.12.36N500d1tol0.001.mat
%                 % this MAT file is to plot figure 3&4 in paper
%             else
%                 warning(['TestcubMCon-gaussian-uniform-Out-17-Aug-2012_13.12.36N500d1tol0.001.mat does not exist. '...
%                     'need to call function RunTestcubMConGaussiand1 to produce the MAT file.'])
%                 [res,test,fun,param] = RunTestcubMConGaussiand1;
            %end
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
    plotTest.logerrlo=-5;
    plotTest.logerrhi=0;
    plotTest.logtimelo=-4;
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
    if any(strcmp('cubSobol',test.whichsample))
        plotTest.err=res.Sobolerr;
        plotTest.time=res.Soboltime;
        plotTest.exit=res.Sobolexit;
        plotTest.name=[plotTest.namepref 'cubSobolErrTime'];
        plotTest.defaultcolor=[1 0 0];
        if any(strcmp('black',plotTest.plotcolor))
            plotTestcubMCblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTestcubMCcolor(plotTest,param)
        end
    end
    
    %% Plot Sobol heavy dutyresults
    if any(strcmp('cubLattice',test.whichsample))
        plotTest.err=res.Latticeerr;
        plotTest.time=res.Latticetime;
        plotTest.exit=res.Latticeexit;
        plotTest.name=[plotTest.namepref 'cubLatticeErrTime'];
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
        %% Plot chebfun and chebfunHeavy results
    if any(strcmp('integral',test.whichsample))
        plotTest.err=res.integralerr;
        plotTest.time=res.integraltime;
        plotTest.name=[plotTest.namepref 'integralErrTime'];
        plotTest.defaultcolor=[1 0 0];
        if any(strcmp('black',plotTest.plotcolor))
            plotTestcubMCblack(plotTest,param)
        end
        if any(strcmp('color',plotTest.plotcolor))
            plotTestcubMCcolor(plotTest,param)
        end
    end
end
end


