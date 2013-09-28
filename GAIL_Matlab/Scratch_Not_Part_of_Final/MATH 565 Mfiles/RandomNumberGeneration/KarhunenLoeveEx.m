%Generating Brownian motion paths
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18)
format compact; clear all; close all
%normcdf = @(x) 0.5 * erfc(-(x) ./ sqrt(2));
%norminv = @(p) -sqrt(2).*erfcinv(2*p);

reset(RandStream.getDefaultStream)

dplot=2.^(0:10)';
dmax=max(dplot);
dvec=(1:dmax)';
zvec=randn(dmax,1);
xplot=0:0.001:1;
nplot=length(xplot);
eigfun=sqrt(2)*sin(((dvec-0.5)*pi)*xplot).*repmat(zvec./(pi*(dvec-0.5)),1,nplot);
brown=cumsum(eigfun,1);
maxbrown=max(max(abs(brown)))*1.1;
for i=1:length(dplot);
    d=dplot(i);
    figure; 
    plot(xplot,brown(d,:),'b-','linewidth',2);
    axis([0 1 -maxbrown maxbrown]);
    title(['{\it d = ' int2str(d) '}'])
    eval(['print -depsc KHd-' int2str(i) '.eps'])
end


