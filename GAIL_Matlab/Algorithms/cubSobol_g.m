function [q,err,time,n,overbudget]=cubSobol(f,d,tol,density)
% [q,err,time,n,overbudget]=CUBSOBOL(f,d,tol) estimates the integral of f 
% over the d-dimensional region using Sobol' cubature

tic
%% Initialize parameters
if nargin < 4
   density='uniform';
   if nargin < 3
      tol=1e-4;
      if nargin < 2
         d=1;
      end
   end
end
if strcmp(density,'normal')
   f=@(x) f(norminv(x));
end
mmax=22; %maximum number of points is 2^mmax
mmin=10; %initial number of points is 2^mmin
mlag=4; %distance between coefficients summed and those computed
fudge=3; %inflation factor
sobstr=sobolset(d); %generate a Sobol' sequence
sobstr=scramble(sobstr,'MatousekAffineOwen'); %scramble it
Stilde=zeros(mmax-mmin+1,1); %initialize
errest=zeros(mmax-mmin+1,1);
appxinteg=zeros(mmax-mmin+1,1);
overbudget=true;

%% Initial points and FWT
n=2^mmin;
xpts=sobstr(1:n,1:d); n0=n;
y=f(xpts);
yval=y;

%% Compute initial FWT
for l=0:mmin-1
   nl=2^l;
   nmminlm1=2^(mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+oddval)/2;
   y(~ptind)=(evenval-oddval)/2;
end
%y now contains the FWT coefficients

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
nllstart=2^(mmin-mlag-1);
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
err=fudge*2^(-mmin)*Stilde(1);
errest(1)=err;
if err <= tol; overbudget=false; time=toc; return, end

%% Loop over m
for m=mmin+1:mmax 
   n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=sobstr(n0+(1:nnext),1:d); 
   n0=n0+nnext;
   ynext=f(xnext);
   yval=[yval; ynext];

   %% Compute initial FWT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+oddval)/2;
      ynext(~ptind)=(evenval-oddval)/2;
   end

   %% Compute FWT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+oddval)/2;
   y(~ptind)=(evenval-oddval)/2;
   %disp('line 100')

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
   nllstart=2^(m-mlag-1);
   meff=m-mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   err=fudge*2^(-m)*Stilde(meff);
   errest(meff)=err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   if err <= tol; overbudget=false; time=toc; return, end

end
time=toc;
