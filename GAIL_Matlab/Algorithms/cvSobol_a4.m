function [q,out_param] = cvSobol_a4(varargin)

tic

%% Check and initialize parameters
[f,g,out_param] = cubSobol_g_param(varargin{:});

if strcmp(out_param.measure,'normal')
   f=@(x) f(gail.stdnorminv(x));
   g=@(x) g(gail.stdnorminv(x));
end

%% Main algorithm
r_lag=out_param.r; %distance between coefficients summed and those computed
l_star=out_param.mmin-r_lag;
sobstr=sobolset(out_param.d); %generate a Sobol' sequence
sobstr=scramble(sobstr,'MatousekAffineOwen'); %scramble it
Stilde=zeros(out_param.mmax-out_param.mmin+1,1); %initialize sum of DFWT terms
StildeNC=zeros(out_param.mmax-out_param.mmin+1,r_lag); %initialize various sums of DFWT terms for necessary conditions
cond1=(1+out_param.fudge(r_lag))*(1+2*out_param.fudge(r_lag-(1:r_lag)))./(1+out_param.fudge(r_lag-(1:r_lag))); % Factors for the necessary conditions
cond2=(1+out_param.fudge(r_lag-(1:r_lag)))*(1+2*out_param.fudge(r_lag))/(1+out_param.fudge(r_lag)); % Factors for the necessary conditions
errest=zeros(out_param.mmax-out_param.mmin+1,1); %initialize error estimates
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1); %initialize approximations to integral
exit_len = 2;
out_param.exit=false(1,exit_len); %we start the algorithm with all warning flags down

%% Initial points and FWT
out_param.n=2^out_param.mmin; %total number of points to start with
n0=out_param.n; %initial number of points
xpts=sobstr(1:n0,1:out_param.d); %grab Sobol' points
% initialize beta
temp = (2^(out_param.mmin-r_lag-1)+1:2^(out_param.mmin));
A=g(xpts);b=f(xpts);beta=L1Reg(A(temp), b(temp));
f1=@(x) f(x) - g(x)*beta;
y=f1(xpts); %evaluate integrand
yval=y;

%% Compute initial FWT
for l=0:out_param.mmin-1
   nl=2^l;
   nmminlm1=2^(out_param.mmin-l-1);
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

%% Create kappanumap implicitly from the data
kappanumap=(1:out_param.n)'; %initialize map
for l=out_param.mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone); %which in the pair are the larger ones
   temp=kappanumap(nl+1+flip); %then flip 
   kappanumap(nl+1+flip)=kappanumap(1+flip); %them
   kappanumap(1+flip)=temp; %around
end

%% Compute Stilde
nllstart=2^(out_param.mmin-r_lag-1);
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
for i = 1:r_lag % Storing the information for the necessary conditions
    nllstart=2*nllstart;
    StildeNC(i,i)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
end
out_param.pred_err=out_param.fudge(out_param.mmin)*Stilde(1);
errest(1)=out_param.pred_err;

deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-errest(1)),...
    out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+errest(1)),out_param.toltype));
deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-errest(1)),...
    out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+errest(1)),out_param.toltype));

is_done = false;
if out_param.pred_err <= deltaplus
   q=q+deltaminus;
   appxinteg(1)=q;
   out_param.time=toc;
   is_done = true;
elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
   out_param.exit(1) = true;
end

%% Loop over m
for m=out_param.mmin+1:out_param.mmax
   if is_done,
       break;
   end
   out_param.n=2^m;
   xpts=sobstr(1:2^m, 1:out_param.d);
   temp = (2^(m-r_lag-1)+1:2^m);
   A=g(xpts);b=f(xpts);beta=L1Reg(A(temp), b(temp));
   f1=@(x) f(x) - g(x)*beta;
   y=f1(xpts); yval=y;
   
   %% Compute initial FWT on next points
   for l=0:m-1
      nl=2^l;
      nmminlm1=2^(m-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      evenval=y(ptind);
      oddval=y(~ptind);
      y(ptind)=(evenval+oddval)/2;
      y(~ptind)=(evenval-oddval)/2;
   end


   %% Update kappanumap
   kappanumap=(1:2^m); %initialize map
   for l=m-1:-1:1
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone);
      temp=kappanumap(nl+1+flip);
      kappanumap(nl+1+flip)=kappanumap(1+flip);
      kappanumap(1+flip)=temp;
   end

   %% Compute Stilde
   nllstart=2^(m-r_lag-1);
   meff=m-out_param.mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   for i = 1:r_lag % Storing the information for the necessary conditions
       nllstart=2*nllstart;
       StildeNC(i+meff-1,i)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   end
   % disp((Stilde(meff)*cond1(1:min(meff-1,r_lag)))>=(StildeNC(meff-1,1:min(meff-1,r_lag)))) % Displaying necessary condition 1 results (1 if satisfied)
   % disp((StildeNC(meff-1,1:min(meff-1,r_lag)).*cond2(1:min(meff-1,r_lag)))>=(Stilde(meff)*ones(1,min(meff-1,r_lag)))) % Displaying necessary condition 2 results (1 if satisfied)
   if ~(all((Stilde(meff)*cond1(1:min(meff-1,r_lag)))>=(StildeNC(meff-1,1:min(meff-1,r_lag))))*   ...
    all((StildeNC(meff-1,1:min(meff-1,r_lag)).*cond2(1:min(meff-1,r_lag)))>=(Stilde(meff)*ones(1,min(meff-1,r_lag))))),
        out_param.exit(2) = true;
   end
   out_param.pred_err=out_param.fudge(m)*Stilde(meff);
   errest(meff)=out_param.pred_err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   
    deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-errest(meff)),...
        out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+errest(meff)),out_param.toltype));
    deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-errest(meff)),...
        out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+errest(meff)),out_param.toltype));
   
   if out_param.pred_err <= deltaplus
      q=q+deltaminus;
      appxinteg(meff)=q;
      out_param.time=toc;
      is_done = true;
   elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
      out_param.exit(1) = true;
   end
end

% Alternative
exit_str=2.^(0:exit_len-1).*out_param.exit;
exit_str(out_param.exit==0)=[];
if numel(exit_str)==0;
    out_param.exitflag=0;
else
    out_param.exitflag=exit_str;
end

% Original
% exit_str='';
% if sum(out_param.exit) == 0
%   exit_str = '0';
% else
%   for i=1:exit_len
%     if out_param.exit(i)==1,
%       if i<exit_len 
%         exit_str = strcat(exit_str,{num2str(i)}, {' '});
%       else
%         exit_str = strcat(exit_str,{num2str(i)});
%       end
%     end
%   end
% end
% out_param.exitflag = exit_str;
% out_param = rmfield(out_param,'exit');

out_param.time=toc;
end


%% Parsing for the input of cubSobol_g
function [f, g, out_param] = cubSobol_g_param(varargin)

% Default parameter values
default.d = 1;
default.abstol  = 1e-4;
default.reltol  = 1e-1;
default.measure  = 'uniform';
default.mmin  = 10;
default.mmax  = 24;
default.fudge = @(m) 5*2.^-m;
default.toltype  = 'max';
default.theta  = 1;
default.r=0;
if numel(varargin)<3
    help cubSobol_g
    warning('MATLAB:cubSobol_g:fdnotgiven',...
        'At least, function f and dimension d need to be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    g = @(x) x-0.5;
    out_param.f=f;
    out_param.g=g;
    out_param.d=1;
else
    f = varargin{1};
    g = varargin{2};
    if ~gail.isfcn(f)
        warning('MATLAB:cubSobol_g:fnotfcn',...
            'The given input f was not a function. Example for f(x)=x^2:')
        f = @(x) x.^2;
	g = @(x) x-0.5;
        out_param.f=f;
	out_param.g=g;
        out_param.d=1;
    else
        out_param.f=f;
	out_param.g=g;
        d = varargin{3};
        if ~isnumeric(d) || ~gail.isposint(d) || ~(d<1112)
            warning('MATLAB:cubSobol_g:dnotposint',...
                'The dimension d must be a positive integer less than 1112. Example for f(x)=x^2:')
            f = @(x) x.^2;
            g = @(x) x-0.5;
            out_param.f=f;
	    out_param.g=g;
            out_param.d=1;
        else
            out_param.d=d;
        end
    end
end;

validvarargin=numel(varargin)>3;
if validvarargin
    in3=varargin(4:end);
    for j=1:numel(varargin)-3
    validvarargin=validvarargin && (isnumeric(in3{j}) ...
        || ischar(in3{j}) || isstruct(in3{j}) || gail.isfcn(in3{j}));
    end
    if ~validvarargin
        warning('MATLAB:cubSobol_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in3=varargin{4};
end

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
  f_addParamVal = @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin   
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.measure = default.measure;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
    out_param.toltype = default.toltype;
    out_param.theta = default.theta;
    out_param.r = default.r;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'g',@gail.isfcn);
    addRequired(p,'d',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'r',default.reltol,@isnumeric);
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
        addOptional(p,'toltype',default.toltype,...
            @(x) any(validatestring(x, {'max','comb'})));
        addOptional(p,'theta',default.theta,@isnumeric);
    else
        if isstruct(in3) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        f_addParamVal(p,'r',default.reltol,@isnumeric);
        f_addParamVal(p,'mmin',default.mmin,@isnumeric);
        f_addParamVal(p,'mmax',default.mmax,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@gail.isfcn);
        f_addParamVal(p,'toltype',default.toltype,...
            @(x) any(validatestring(x, {'max','comb'})));
        f_addParamVal(p,'theta',default.theta,@isnumeric);
    end
    parse(p,f,g,d,varargin{4:end})
    out_param = p.Results;
end;

% Force absolute tolerance greater than 0
if (out_param.abstol <= 0 )
    warning('MATLAB:cubSobol_g:abstolnonpos',['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol < 0) || (out_param.reltol > 1)
    warning('MATLAB:cubSobol_g:reltolnonunit',['Relative tolerance should be chosen in [0,1].' ...
            ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Force measure to be uniform or normal only
if ~(strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') )
    warning('MATLAB:cubSobol_g:notmeasure',['The measure can only be uniform or normal.' ...
            ' Using default measure ' num2str(default.measure)])
    out_param.measure = default.measure;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('MATLAB:cubSobol_g:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 54
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning('MATLAB:cubSobol_g:wrongmmax',['The maximum exponent for the budget should be an integer smaller or equal to 53.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('MATLAB:cubSobol_g:fudgenofcn',['The fudge factor should be a positve function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force toltype to be max or comb
if ~(strcmp(out_param.toltype,'max') || strcmp(out_param.toltype,'comb') )
    warning('MATLAB:cubSobol_g:nottoltype',['The error type can only be max or comb.' ...
            ' Using default toltype ' num2str(default.toltype)])
    out_param.toltype = default.toltype;
end

% Force theta to be in [0,1]
if (out_param.theta < 0) || (out_param.theta > 1)
    warning('MATLAB:cubSobol_g:thetanonunit',['Theta should be chosen in [0,1].' ...
            ' Using default theta ' num2str(default.theta)])
    out_param.theta = default.theta;
end

end



