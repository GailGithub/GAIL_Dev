function Latticesuccess = RunTestCubatureonGeoAsianCallLattice
clear all, close all
format compact

fun.funtype='geomean';
fun.S0=100;
fun.K=100;
fun.T=1;
fun.r=0.03;
param.measure='normal';
param.abstol=2e-2;
param.reltol=0; % 0 reltol means all absolut error
param.toltype  = 'max';
param.fudge  = @(m) ones(size(m));
param_indicator=1;

test.nrep=500;
test.howoftenrep=10;
dimchoice=[1 2 4 8 16 32 64]';
ndim=size(dimchoice,1);
test.randch.dim=dimchoice(randi(ndim,test.nrep,1));
sigmin=0.1;
sigmax=0.7;
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
    gail.print_iterations(irep, 'irep', true);
    param.sample='qmc';
    [testfunqmc,fun,param]= ...
          test.randchoicefun(fun,param,test.randch,irep);
    res.dim(irep)=param.dim;
    
    % Evaluate integral using cubLattice
    hyperbox = [-inf(1,param.dim);inf(1,param.dim)];
        [q,out_param]=...
           cubLattice_g(testfunqmc,hyperbox,...
           'measure',param.measure,'abstol',param.abstol,'reltol',param.reltol,'transform','Baker');
        %res.Latticeexit(irep)=cellstr(out_param.exitflag); % we have problems saving this info
        res.LatticeQ(irep)=q;
        res.Latticeerr(irep)=abs(param.exactintegral-q)/gail.tolfun(param.abstol,...
    param.reltol,0,param.exactintegral,param.toltype)*param_indicator;
        res.Latticetime(irep)=out_param.time;
        res.Latticeneval(irep)=out_param.n;
end
% timestamp=datestr(now,'yyyy-mm-dd-HH-MM');
% save(['TestCubature-' fun.funtype '-' param.measure '-' timestamp ...
%    '-N-' int2str(test.nrep) '-d-' int2str(param.dim)  ...
%     '-tol-' num2str(param.abstol) '.mat'])


%% Display test results
%Display TestDiffSettings results
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','tex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markersset(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
plotTest.logerrlo=-6;
plotTest.logerrhi=1;

param.tol=param_indicator;
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
    plotTestColor(plotTest,param)
    end
    Latticesuccess=mean(res.Latticeerr<=param_indicator);
end
gail.save_mat('Paper_cubLattice_g', 'Paper_cubLattice_g_TestGeoAsianCall', true, Latticesuccess,dimchoice,...
        fun,irep,res,test,testfunqmc);    
end

%%%%% Above is all execution. Below is just function definitions

%% Defining function plotTestColor
function plotTestColor(plotTest,param)
[~,~,MATLABVERSION] = GAILstart(false);
if usejava('jvm') || MATLABVERSION <= 7.12
figure
ntot=length(plotTest.err);
colorscatter=repmat(plotTest.defaultcolor,ntot,1);
if isfield(plotTest,'kurtvec') && isfield(plotTest,'kurtmax')
    smallkurt=plotTest.kurtvec<=plotTest.kurtmax;
    colorscatter(smallkurt,:)=repmat([0 0 1],sum(smallkurt),1);
end
if isfield(plotTest,'exit')
    bigsample=plotTest.exit==1;
    colorscatter(bigsample,:)=repmat([0.65 0.35 0],sum(bigsample),1);
end
erraug=[plotTest.errlowlimit; sort(plotTest.err); plotTest.errhilimit];
timeaug=[plotTest.timelowlimit; sort(plotTest.time); plotTest.timehilimit];
probaug=[0; ((0:ntot-1)'+1/2)/ntot; 1];
scatter(plotTest.err,plotTest.time,plotTest.ptsize,colorscatter,'.');
ax1 = gca;
set(ax1,'XColor','k','YColor','k',...
    'XLim',[plotTest.errlowlimit plotTest.errhilimit],...
    'YLim',[plotTest.timelowlimit plotTest.timehilimit],...
    'LineWidth',plotTest.linewidth, ...
    'XScale','log','Yscale','log', ...
    'XTick',10.^(plotTest.logerrlo:plotTest.logerrhi))
    %'XGrid','on','YGrid','on','GridLineStyle','--',...
line([param.tol param.tol],[plotTest.timelowlimit plotTest.timehilimit],...
    'color','k','linestyle','--','linewidth',plotTest.linewidth)
wherexsuccess=-0.02+0.75*(log10(param.tol)-plotTest.logerrlo)/...
    (plotTest.logerrhi-plotTest.logerrlo);
wherexfailure=0.28+0.75*(log10(param.tol)-plotTest.logerrlo)/...
    (plotTest.logerrhi-plotTest.logerrlo);
whereysuccess=0.2;
whereyfailure=0.85;
annotation('textarrow',[wherexsuccess,0.15],whereysuccess*[1 1],...
    'String','Success','FontWeight','bold',...
    'linewidth',2);
annotation('textarrow',[wherexfailure,0.88],whereyfailure*[1 1],...
    'String','Failure','FontWeight','bold',...
    'linewidth',2);
xlabel('Error')
ylabel('Time (seconds)')
axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor',[0 0.7 0],...
           'XLim',[plotTest.errlowlimit plotTest.errhilimit],...
           'Xscale','log','Xtick',[],...
           'YLim',[0 1],'Linewidth',plotTest.linewidth);
%ylabel('Probability')
line(erraug,probaug,'color',[0 0.7 0],'linewidth',plotTest.linewidth)
axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','m','YColor','k',...
           'YLim',[plotTest.timelowlimit plotTest.timehilimit],...
           'Yscale','log','Ytick',[],...
           'XLim',[0 1],'Linewidth',plotTest.linewidth);
%xlabel('Probability')
line(probaug,timeaug,'color','m','linewidth',plotTest.linewidth)
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_TestGeoAsianCall');
%print('-depsc',[plotTest.name '.eps'])
% print('-depsc', ['./Results/' plotTest.name '.eps'])
close all
end
end

%% Random choice GeoCall
function [testfun,fun,out_param]=randchoiceGeoCall(fun,in_param,rchparam,irep)
out_param=in_param;
out_param.dim=rchparam.dim(irep);
out_param.interval=[-Inf(1,out_param.dim); Inf(1,out_param.dim)];
%keyboard
[testfun,out_param]=geomMeanAsianCall(fun,out_param);
end

%% GeoMean Asian Call
function [testfun,param]=geomMeanAsianCall(fun,param)
%   This function chooses and sets up a test function from the parameters
%      input by the user and contained in the structures
%      fun and param
%   fun.funtype         = type of test function
%   param.interval      = domain of test function
%   param.dim           = dimension of the domain
%   param.measure       = probability density function for integration
%   param.exactintegral = exact value of the integral (scalar)
%   fun.shape           = shape parameter (1 x param.dim)
%   fun.scale           = scale parameter (1 x param.dim)
%   fun.addc            = additive constant (1 x param.dim)
%   fun.overaddc        = overall additive constant (scalar)
%   fun.overmultc       = overall multiplicative constant (scalar)

if nargin < 2; %give the basic default parameters 
    param.interval=[0;1]; %default integration interval
    if nargin < 1; fun.funtype='exp'; end %exponential test function
end
if ~isfield(param,'interval'); param.interval=[0;1]; end %default interval

%[~,param]=cubMCparam([],param,'fun'); %check validity of some parameters
%keyboard

if strcmp(fun.funtype,'geomean') %geometric mean Asian call test function
    [testfun,param]=makeGeometricMeanTestFun(fun,param);
else
    error('Function type not recognized')
end
end

%% Geometric Mean Call Integrand
function [testfun,param]=makeGeometricMeanTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'S0','K','T','d','r','sigma'}, ...
    {[1 1],[1 1],[1 1],[1 1],[1 1],[1 1]}, ...
    {100,100,1,1,0.03,0.5});
param.fun=fun; %copy of the function parameters

%if false %Force bb
if strcmp(param.sample,'iid') %Time discretization calculations
   rT=fun.r*fun.T;
   rTfac=((fun.r-fun.sigma^2/2)*fun.T/param.dim)*(1:param.dim);
   Delta=fun.sigma*sqrt(fun.T/param.dim);
   %keyboard
   testfun=@(x) geomeancalltd(x,fun.S0,fun.K,rT,param.dim,rTfac,Delta);
else %Brownian bridge calculations for qMC
   log2d=log2(param.dim);
   rT=fun.r*fun.T;
   rTfac=((fun.r-fun.sigma^2/2)*fun.T/param.dim)*(1:param.dim);
   dov2l=param.dim;
   ramp=zeros(log2d+1,param.dim);
   ramp(1,:)=(fun.sigma*sqrt(fun.T)/dov2l)*(1:dov2l);
   twolmin1=1;
   Deltal=fun.T/4;
   dov2lmin1=zeros(log2d,1);
   for l=1:log2d
       dov2lmin1(l)=dov2l;
       dov2l=dov2l/2;
       ramppc=(fun.sigma*sqrt(Deltal)/dov2l)*[1:dov2l dov2l-1:-1:0];
       ramp(l+1,:)=repmat(ramppc,1,twolmin1);
       Deltal=Deltal/2;
       twolmin1=twolmin1*2;
   end
   testfun=@(x) geomeancallbb(x,fun.S0,fun.K,rT,param.dim,log2d,rTfac,ramp,dov2lmin1);
end

%% Compute exact integral of this function
rTmod=(fun.r-fun.sigma^2/2)*fun.T/2;
sigrootTmod=fun.sigma*sqrt(fun.T*(1-1/(2*param.dim)^2)/3);
x0=(log(fun.K/fun.S0)-rTmod)/sigrootTmod;
param.exactintegral=exp(-fun.r*fun.T)*(fun.S0*exp(rTmod+sigrootTmod^2/2)...
    *normcdf(sigrootTmod-x0) - fun.K*normcdf(-x0));
end

function payoff=geomeancalltd(x,S0,K,rT,d,rTfac,Delta)
    n=size(x,1);
    %Create stock paths
    %keyboard
    Smatrix=repmat(rTfac,n,1)+cumsum(Delta*x,2);
    %Compute payoff
    loggmean=sum([Smatrix(:,1:d-1) Smatrix(:,d)/2],2)/d;
    %keyboard
    payoff=max(S0*exp(loggmean)-K,zeros(n,1))*exp(-rT);
end 

function payoff=geomeancallbb(x,S0,K,rT,d,log2d,rTfac,ramp,dov2lmin1)
    n=size(x,1);
    %Create stock paths
    Smatrix=repmat(rTfac,n,1);
    Smatrix=Smatrix+x(:,1)*ramp(1,:);
    twolmin1=1;
    for l=1:log2d
        jvec=twolmin1+(1:twolmin1);
        %keyboard
        weight=reshape(repmat(x(:,jvec),dov2lmin1(l),1),n,d);
        %keyboard
        Smatrix=Smatrix+weight.*repmat(ramp(l+1,:),n,1);
        twolmin1=twolmin1*2;
    end
    %Compute payoff
    loggmean=sum([Smatrix(:,1:d-1) Smatrix(:,d)/2],2)/d;
    %keyboard
    payoff=max(S0*exp(loggmean)-K,zeros(n,1))*exp(-rT);
end
