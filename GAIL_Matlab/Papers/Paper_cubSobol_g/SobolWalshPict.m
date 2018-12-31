function SobolWalshPict
%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
clear all %clear all variables
close all %close all figures
%Some defaults to make the plots easier to view
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultLineMarkerSize',40) %large dots

% %% Plot Walsh functions
% mmax=3;
% kmax=2^mmax;
% invkmax=1/kmax;
% hplot=0.002;
% xplot=bsxfun(@plus,[(0:hplot:invkmax-hplot)'; invkmax-100*eps],...
%    (0:invkmax:1-invkmax));
% xplot=xplot(:);
% nplot=length(xplot);
% walshplot=ones(nplot,kmax);
% for m=1:mmax;
%    n=2^m;
%    whichchange=mod(floor((0:kmax-1)*(2/n)),2)==1;
%    %keyboard
%    walshplot(:,whichchange)=walshplot(:,whichchange).*...
%       repmat((-1).^floor(n*xplot),1,kmax/2);
% end
% for k=0:kmax-1;
%    figure
%    plot(xplot,walshplot(:,k+1),'b-')
%    xlabel('$x$')
%    ylabel(['walsh($' int2str(k) ',x)$'])
%    axis([0 1 -1.2 1.2])
%    set(gca,'xtick',0:0.25:1,'ytick',-1:0.5:1)
%    % eval(['print -depsc walshk' int2str(k) 'fun.eps'])
% end

[~,~,MATLABVERSION] = GAILstart(false);
if usejava('jvm') || MATLABVERSION <= 7.12

%% Plot Sobol points
mvec=(8:8); %powers of 2
nm=length(mvec); %number of plots
nvec=2.^mvec; %numbers of points to plot
nmax=nvec(end); %maximum number of points
rng(61679); %to get the same answer each time
d=2;
pscsob=sobolset(d);
xsob=net(pscsob,nmax);
for j=1:nm
   n=nvec(j);
   figure
   plot(xsob(1:n,1),xsob(1:n,2),'k.')
      %xsob([2 6 5],1),xsob([2 6 5],2),'r.')
%    xlabel('$x_1$')
%    ylabel('$x_2$')
   axis equal
   axis([0 1 0 1])
   set(gca,'xtick',0:0.25:1,'ytick',0:0.25:1)
%    text(xsob(2,1)+0.03,xsob(2,2),'$\mathbf{z}_1$','BackgroundColor','white')
%    text(xsob(6,1)+0.03,xsob(6,2),'$\mathbf{z}_5$','BackgroundColor','white')
%    text(xsob(5,1)+0.03,xsob(5,2),'$\mathbf{z}_4$','BackgroundColor','white')
   %eval(['print -depsc sob' int2str(n) 'pts.eps'])
   gail.save_eps('Paper_cubSobol_g', 'Paper_cubSobol_g_256SobolPoints');
end

%% Plot scrambled Sobol points
mvec=(8:8); %powers of 2
nm=length(mvec); %number of plots
nvec=2.^mvec; %numbers of points to plot
nmax=nvec(end); %maximum number of points
rng(61679); %to get the same answer each time
d=2;
pscsob=scramble(sobolset(d),'MatousekAffineOwen');
xsob=net(pscsob,nmax);
for j=1:nm
   n=nvec(j);
   figure
   plot(xsob(1:n,1),xsob(1:n,2),'k.')
      %xsob([2 6 5],1),xsob([2 6 5],2),'r.')
%    xlabel('$x_1$')
%    ylabel('$x_2$')
   axis equal
   axis([0 1 0 1])
   set(gca,'xtick',0:0.25:1,'ytick',0:0.25:1)
%    text(xsob(2,1)+0.03,xsob(2,2),'$\mathbf{z}_1$','BackgroundColor','white')
%    text(xsob(6,1)+0.03,xsob(6,2),'$\mathbf{z}_5$','BackgroundColor','white')
%    text(xsob(5,1)+0.03,xsob(5,2),'$\mathbf{z}_4$','BackgroundColor','white')
gail.save_eps('Paper_cubSobol_g', 'Paper_cubSobol_g_256ScrambledShiftedSobolPoints');
   %eval(['print -depsc scrsob' int2str(n) 'pts.eps'])
end
close all
end

%  Construct an MException object to represent the error.
msgID = 'MYFUN:DemoExample';
msg = 'Demonstrate throw exception.';
baseException = MException(msgID,msg);

% Throw the exception to stop execution and display an error
% message.
%throw(baseException)

end
