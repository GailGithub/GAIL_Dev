%Automatic cubature with Sobol sequences

%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',30,'defaulttextfontsize',30) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %larger markers
tic

%% Initialize parameters
tol=1e-4;
%testfun=@(x) x; exactinteg=1/2; d=1; %test function
%testfun=@(x) x.^2; exactinteg=1/3; d=1; %test function
%a=20; testfun=@(x) sin(a*x); exactinteg=(1-cos(a))/a; d=1; %test function
%testfun=@(x) x(:,1).*x(:,2); exactinteg=1/4; d=2; %test function
%testfun=@(x) sin(x(:,1)).*x(:,2)+exp(x(:,1)); exactinteg=(1-cos(1))/2 + (exp(1)-1); d=2; %test function
a=3; d=5; testfun=@(x) exp(a*sum(x,2))./(((exp(a)-1)/a).^d); exactinteg=1; %test function

[q,errest,time]=cubSobol(testfun,d,tol)
trueerr=abs(exactinteg-q)
break

%% Plot results
figure
mfinal=m;
minexp=floor(mmin*log10(2));
maxexp=ceil(mfinal*log10(2));
h=loglog(2.^(mmin:mfinal),trueerr(1:meff),'b.',2.^(mmin:mfinal),errest(1:meff),'rs');
set(h(2),'MarkerFaceColor','r','MarkerSize',10);
set(gca,'Xtick',10.^(minexp:maxexp))
ymin=10.^floor(log10(max(1e-15, min([trueerr; errest]))));
axis([10.^[minexp maxexp] ymin 1])



