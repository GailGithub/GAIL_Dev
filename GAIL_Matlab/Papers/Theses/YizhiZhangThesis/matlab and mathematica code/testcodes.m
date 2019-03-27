t=0.1;
delta=0.2; 
f = @(x) (1/delta.^4).*((x-t).^3/6.*(x>=t).*(x<t+delta)...
        +(-3.*(x-t).^3+12.*delta.*(x-t).^2-12.*delta.^2.*(x-t)+4.*delta.^3)/6.*(x>=t+delta).*(x<t+2*delta)...
        +(3.*(x-t).^3-24.*delta.*(x-t).^2+60.*delta.^2.*(x-t)-44.*delta.^3)/6.*(x>=t+2*delta).*(x<t+3*delta)...
        +(t+4.*delta-x).^3/6.*(x>=t+3*delta).*(x<=t+4*delta));

% f = @(x) exp(-x.^8);
    
avec=linspace(0,10,1000);
nrep=length(avec);
ntrapvec=zeros(nrep,1);
nsimpvec=ntrapvec;
tolvec=ntrapvec;
errtrapvec=ntrapvec;
errsimpvec=ntrapvec;

for i=1:nrep
    tolvec(i)=10.^(-avec(i));
    [qtrap,out_paramt]=integral_t(f,'abstol',tolvec(i),'nmax',1e7,'hcut',0.1);
    ntrapvec(i)=out_paramt.npoints;
    errtrapvec(i)=abs(1-qtrap);
    [qsimp,out_params]=integral_s(f,'abstol',tolvec(i),'nmax',1e7,'hcut',0.1);
    nsimpvec(i)=out_params.npoints;
    errsimpvec(i)=abs(1-qsimp);
end
%     disp(epsilonvec);
slopehalfvec=tolvec.^(-0.5);
slopequarvec=tolvec.^(-0.25);

figure
loglog(tolvec,ntrapvec)
hold on
loglog(tolvec,nsimpvec,'r')
hold on
loglog(tolvec,slopehalfvec,'g')
hold on
loglog(tolvec,slopequarvec,'k')
hold off
legend('trapezoidal','Simpson','slope -1/2','slope -1/4','location','northeast')
title('convergence rate for two algorithms')
xlabel('absolute tolerance')
ylabel('number of points')