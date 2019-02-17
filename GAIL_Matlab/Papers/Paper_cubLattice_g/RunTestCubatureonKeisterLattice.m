function Latticesuccess = RunTestCubatureonKeisterLattice
clearvars, close all
format compact

fun.funtype='Keister';
param.measure='uniform';
param.abstol=2e-2; % 0 abstol means all relative error
param.reltol=0; % 0 reltol means all absolute error
param.toltype  = 'max';
param_indicator=1;

test.nrep=500;
test.howoftenrep=10;
test.randch.a=sqrt(2);
test.randch.dim=floor(20.^rand(test.nrep,1)); %random dimensions to 20
test.randchoicefun=@randchoiceKeister;
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
    hyperbox = [zeros(1,param.dim);ones(1,param.dim)];
        [q,out_param]=...
           cubLattice_g(testfunqmc,hyperbox,...
           'measure',param.measure,'abstol',param.abstol,'reltol',param.reltol,'transform','Baker');
        %res.Latticeexit(irep)=cellstr(out_param.exitflag); % we have problems saving this info
        res.LatticeQ(irep)=q;
        res.Latticeerr(irep)=abs(param.exactintegral-q)/gail.tolfun(param.abstol,...
    param.reltol,0,param.exactintegral,param.toltype)*param_indicator;
        res.Latticetime(irep)=out_param.time;
        res.Latticeneval(irep)=out_param.n;

% timestamp=datestr(now,'yyyy-mm-dd-HH-MM');
% save(['TestCubature-' fun.funtype '-' param.measure '-' timestamp ...
%    '-N-' int2str(test.nrep) '-d-' int2str(param.dim)  ...
%     '-tol-' num2str(param.abstol) '.mat'])
end

%% Display test results
%Display TestMCDiffSettings results
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','tex') %tex axis labels
set(0,'defaultLineMarkerSize',40) %larger markersset(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
plotTest.logerrlo=-4;
plotTest.logerrhi=2;

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
gail.save_mat('Paper_cubLattice_g', 'Paper_cubLattice_g_TestKeister',true, Latticesuccess,...
        fun,irep,res,test);
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
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_TestKeister');
%print('-depsc',[plotTest.name '.eps'])
% print('-depsc', ['./Results/' plotTest.name '.eps'])
close all
end
end

%% Random choice Keister
function [testfun,fun,param]=randchoiceKeister(fun,param,rchparam,irep)
%Choose Keister function
param.a=rchparam.a;
param.dim=rchparam.dim(irep);
param.interval=[zeros(1,param.dim); ones(1,param.dim)];
param.exactintegral = Keistertrue(param.dim);
testfun = @(x) Keisterfun(x,param.dim,param.a)./param.exactintegral; 
   %normalize to unity
param.exactintegral=1;
%keyboard
end

function f=Keisterfun(x,d,a) %a must be bigger than sqrt(2)
sumsq=sum(norminv(x(:,1:d)).^2,2);
f=((sqrt(2*pi)/a).^d).*exp(-(1./(a.^2)-1./2)*sumsq).*cos(sqrt(sumsq)/a);
%keyboard
end

function I = Keistertrue(d)
%KEISTERTRUE computes the true value of the Keister integral in dimension d
%  accuracy might degrade as d increases due to round-off error
cosinteg=zeros(1,d);
cosinteg(1)=sqrt(pi)/(2*exp(1/4));
sininteg=zeros(1,d);
%sininteg(1)=integral(@(x) exp(-x.*x).*sin(x),0,inf);
sininteg(1)=4.244363835020225e-01;
cosinteg(2)=(1-sininteg(1))/2;
sininteg(2)=cosinteg(1)/2;
for j=3:d
   cosinteg(j)=((j-2)*cosinteg(j-2)-sininteg(j-1))/2;
   sininteg(j)=((j-2)*sininteg(j-2)+cosinteg(j-1))/2;
end
I=(2*(pi.^(d/2))/gamma(d/2))*cosinteg(d);
end


