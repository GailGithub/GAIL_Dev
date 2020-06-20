%This script will use the beginfool() function to find points called for by
%the function of your choosing (quad, quadgk, chebfun, etc.)
%It then fits peaks in between the points it collects
%Just change this script, and you will get a test of beginfool AND peakyfunction

%This script also gives the ratio of the max derivative values. 
%The last addition was finding the "real" integral so we could get an error

%% Garbage collection
function FoolAutomaticAlgorithms()
clearvars, close all, format long, format compact
set(0,'defaulttextfontsize',24,'defaultaxesfontsize',24)
hold all

info.filename='scriptValues.mat'; %This file name will be passed throughout the script
xsample=[];
saveMCQMC2012peakyfundir(info.filename,xsample)
%% ---------------Function and Bounds----------
%Pick the function you want to put peaks into.
info.degree=2;
info.coefficient=0;
info.RegFunc=... %slowly varying function
    @(x) info.coefficient*x.^info.degree;
info.RegFuncprime=... %and its derivative
    @(x) info.degree.*info.coefficient*x.^(info.degree-1); %slowly varying
info.RegFuncdubprime=... %and its second derivative
    @(x) info.degree.*(info.degree-1).*info.coefficient*x.^(info.degree-2) 
info.lower=0; %left endpoint of function domain
info.upper=1; %right endpoint of function domain
info.p=5; %peak smoothness
info.c=1;
info.sign=1; %sign of the peaks

%% Ways to call a function
fname='chebint' %insert the automatic algorithm that you want to call
switch fname %give the calling sequence
    case 'quadgk'
        callautoalg = @(fun,lower,upper) quadgk(fun,lower,upper,'AbsTol',1e-14);
    case 'quad'
        callautoalg = @(fun,lower,upper) quad(fun,lower,upper,1e-14);
    case 'chebint'        
        if exist('chebfun','file') 
          callautoalg = @(fun,lower,upper) sum(chebfun(fun,[lower upper]));
        else 
           warning('Please install the chebfun package to plot this figure.')
        end
    case 'fminbnd'
        callautoalg = @(fun,lower,upper) fminbnd(fun,lower,upper);
    case '???' %put the NAG algorithm here
        callautoalg = @(fun,lower,upper) fminbnd(fun,lower,upper);
    case 'integral'
        callautoalg = @(fun,lower,upper) integral(fun,lower,upper,'AbsTol',1e-14);
end
 

%% Call the automatic algorithm to fool
tic
if exist('chebfun','file') 
  Original=callautoalg(@(x) snooper(x,info),info.lower,info.upper)
end
%% Construct the peaks with zeros at the points called by the
%  auotmatic algorithm
peaks=@(x) peakyfunction(x,info);

%% -----------------Fooling the Automatic Algorithm----------
%The following function needs to be the one you are tricking
if exist('chebfun','file') 
  inaccurate=callautoalg(peaks,info.lower,info.upper);
end

%if 'inaccurate' = 'Original' then peakyfunction successfully broke quad, or
%quadgk, etc.
toc

%% ----------------Plotting the Peaky Function-------
x=info.lower:.002:info.upper; %The domain of the peaky function
[yplot, primeplot, dubplot, info]=peaks(x);

%Here the peaks are plotted over the original function
subplot(2,3,1), plot(x,yplot,'b-',info.sortedX,info.RegFunc(info.sortedX),'ro'),title('Peaky Function Overlaid'),...
    %legend('Peaky Function','Original','Location','EastOutside')

%plots of Derivatives
 subplot(2,3,3), plot(x,primeplot),title('First Derivative') 
 subplot(2,3,5), plot(x,dubplot),title('Second Derivative')  

%% --------------------Working with the Derivatives-------------
info.secondzeros=[info.sortedX(2:end)'; info.sortedX(1:end-1)'; (.5*((info.sortedX(2:end)...
    -info.sortedX(1:end-1))/sqrt(2*info.p-1)+info.sortedX(1:end-1)+...
    info.sortedX(2:end)))'; (.5*((-info.sortedX(2:end)+info.sortedX(1:end-1))/...
    sqrt(2*info.p-1)+info.sortedX(1:end-1)+info.sortedX(2:end)))'];
[~,yp]=peaks(info.secondzeros);
primemaxvec=max(abs(yp));
info.thirdzeros=[info.sortedX(2:end)'; info.sortedX(1:end-1)'; ((info.sortedX(1:end-1)+...
    info.sortedX(2:end))/2)'; (.5*(sqrt(3)*(info.sortedX(2:end)-info.sortedX(1:end-1))/...
    sqrt(2*info.p-1)+info.sortedX(1:end-1)+info.sortedX(2:end)))'; (.5*(sqrt(3)*(-info.sortedX(2:end)...
    +info.sortedX(1:end-1))/sqrt(2*info.p-1)+info.sortedX(1:end-1)+info.sortedX(2:end)))'];
[~,~,ydub]=peaks(info.thirdzeros);
dubmaxvec=max(abs(ydub));
primemax=max(primemaxvec);
dubmax=max(dubmaxvec);

ratio=dubmax/primemax

%% -----------------------Real Integral(for p=integer?)----------------------
upperbnd=info.sortedX(2:end);
lowerbnd=info.sortedX(1:end-1);
withbumps=(upperbnd-lowerbnd)/2*(sqrt(pi)*gamma(info.p+1)/gamma(1.5+info.p));
[~,~,MATLABVERSION] = GAILstart(false);
if exist('chebfun','file') && (usejava('jvm') || MATLABVERSION <= 7.12)
    failintegral=inaccurate
    realintegral=info.c*(sum(withbumps))+info.coefficient/(info.degree+1)*...
        (info.upper^(info.degree+1)-info.lower^(info.degree+1))
    %relative error
    error=abs((realintegral-failintegral)/realintegral)
    hold off
    
    figure
    if strcmp(fname,'quadgk') || strcmp(fname,'integral')%too many peaks
        x=info.lower:.00001:info.upper*0.01;
        yplot=peaks(x);
        axisright=info.upper*0.01;
    else
        axisright=info.upper;
    end
    %This one is for slides
    plot(x,yplot/realintegral,'b-',...
        info.sortedX,peaks(info.sortedX)/realintegral,'r.',...
        'linewidth',3,'markersize',30)
    axis([info.lower axisright 0 3])
    gail.save_eps('MCQMC2012PaperOutput/Results',[fname 'color']);
    %eval(['print -depsc Fool' fname 'color.eps'])
    %This one is for printing in articles
    plot(x,yplot/realintegral,'k-',...
        info.sortedX,peaks(info.sortedX)/realintegral,'k.',....
        'linewidth',3,'markersize',30)
    axis([info.lower axisright 0 3])
    gail.save_eps('MCQMC2012PaperOutput/Results',[fname 'bw']);
    %eval(['print -depsc Fool' fname 'bw.eps'])
end
end


