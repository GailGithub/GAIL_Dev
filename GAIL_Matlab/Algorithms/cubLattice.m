function [q,err,time,n]=cubLattice(f,d,tol)
% [q,err,time]=cubLattice(f,d,tol) estimates the integral of f over the
% d-dimensional unit cube using rank-1 Lattice rules
% f is the function predefined.
% d is the dimension of the function.
% tol is our given error tolerance.
% Example:
%   f=@(x) x(:,1).*x(:,2);
%   d=2;
%   tol=10^-5;
%   trueint=1/4;
%   error=abs(trueint-cubLattice(f,d,tol));

tic
%% Initialize parameters
if nargin < 3
   tol=1e-4;
   if nargin < 2
      d=1;
   end
end
mmax=17; %maximum budget of points is 2^mmax. Check with the lattice generator.
mmin=10; %initial number of points is 2^mmin
mlag=4;
fudge=3;
Stilde=zeros(mmax-mmin+1,1);
errest=zeros(mmax-mmin+1,1);
appxinteg=zeros(mmax-mmin+1,1);
shift=rand; % Random shift to our lattice rule in order to avoid discontinuities.

%% Applying the Baker's transform to make any input function periodic
f=@(x) f(1-2*abs(x-1/2));

%% Initial points and FWT
n=2^mmin;
xpts=mod(lattice_gen(1,n,d)+shift,1);
n0=n;
y=f(xpts);
yval=y;

%% Compute initial FFT
for l=0:mmin-1
   nl=2^l;
   nmminlm1=2^(mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;
end

%% Approximate integral
q=mean(yval);
appxinteg(1)=q;

%% Create kappanumap
kappanumap=(1:n)'; %initialize map
for l=mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone); %
   temp=kappanumap(nl+1+flip);
   kappanumap(nl+1+flip)=kappanumap(1+flip);
   kappanumap(1+flip)=temp;
end

%% Compute Stilde
nllstart=int64(2^(mmin-mlag-1));
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
err=fudge*2^(-mmin)*Stilde(1);
errest(1)=err;
if err <= tol; time=toc; return, end

timeT=0;
%% Loop over m
for m=mmin+1:mmax 
   n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=mod(lattice_gen(nnext+1,2*nnext,d)+shift,1);
   n0=n0+nnext;
   ynext=f(xnext);
   yval=[yval; ynext];

   %% Compute initial FFT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
      coefv=repmat(coef,nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+coefv.*oddval)/2;
      ynext(~ptind)=(evenval-coefv.*oddval)/2;
   end

   %% Compute FFT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;

   %% Update kappanumap
   kappanumap=[kappanumap; (nnext+1:n)']; %initialize map
   for l=m-1:-1:1
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone); %
      temp=kappanumap(nl+1+flip);
      kappanumap(nl+1+flip)=kappanumap(1+flip);
      kappanumap(1+flip)=temp;
   end

   %% Compute Stilde
   nllstart=int64(2^(m-mlag-1));
   meff=m-mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   err=fudge*2^(-m)*Stilde(meff);
   errest(meff)=err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   if err <= tol; time=toc; break
    elseif m==mmax; warning('The maximum budget is attained without reaching the tolerance.');
   end

end
time=toc;

