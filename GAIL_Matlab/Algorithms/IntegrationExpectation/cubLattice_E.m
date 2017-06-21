function [q,out_param,y,kappanumap] = cubLattice_g(varargin)

t_start = tic;
%% Initial important cone factors and Check-initialize parameters
r_lag = 4; %distance between coefficients summed and those computed
[f,hyperbox,out_param, cv, FuncCount, CVCount] = cubLattice_g_param(r_lag,varargin{:});

%------------------------------------------------------------------------------
% TRANSFORMATION
%changing the integrand and the hyperbox when measure is uniform ball or
%sphere by applying the appropriate transformation
if strcmpi(out_param.measure,'uniform ball') || strcmpi(out_param.measure,'uniform sphere')% using uniformly distributed samples on a ball or sphere
    if strcmp(out_param.measure,'uniform sphere') && out_param.transf == 1 %box-to-sphere transformation
        out_param.d = out_param.d + 1; % changing out_param.d to the dimension of the sphere
        out_param.shift = [out_param.shift rand];
    end
    
    if strcmpi(out_param.measure,'uniform ball')% using the formula of the volume of a ball
        volume = ((2.0*pi^(out_param.d/2.0))/(out_param.d*gamma(out_param.d/2.0)))*out_param.radius^out_param.d; %volume of a d-dimentional ball
    else % using the formula of the volume of a sphere
        volume = ((2.0*pi^(out_param.d/2.0))/(gamma(out_param.d/2.0)))*out_param.radius^(out_param.d - 1); %volume of a d-dimentional sphere
    end
    
    if out_param.transf == 1 % box-to-ball or box-to-sphere transformation should be used
        if out_param.d == 1 % It is not necessary to multiply the function f by the volume, since no transformation is being made
            hyperbox = [hyperbox - out_param.radius; hyperbox + out_param.radius];% for one dimension, the ball is actually an interval
            out_param.measure = 'uniform';% then a uniform distribution on a box can be used
        else
            if strcmpi(out_param.measure,'uniform ball') % box-to-ball transformation
                f = @(t) f(gail.domain_balls_spheres.ball_psi_1(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
            else %  % box-to-sphere transformation
                f = @(t) f(gail.domain_balls_spheres.sphere_psi_1(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
                out_param.d = out_param.d - 1;% the box-to-sphere transformation takes points from a (d-1)-dimensional box to a d-dimensional sphere
                out_param.shift = out_param.shift(1:end-1);
            end
            hyperbox = [zeros(1, out_param.d); ones(1, out_param.d)];% the hyperbox must be the domain of the transformation, which is a unit box
            out_param.measure = 'uniform';% then a uniform distribution on a box can be used
        end
    else % normal-to-ball or normal-to-sphere transformation should be used
        if strcmpi(out_param.measure,'uniform ball') % normal-to-ball transformation
            f = @(t) f(gail.domain_balls_spheres.ball_psi_2(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
        else % normal-to-sphere transformation
            f = @(t) f(gail.domain_balls_spheres.sphere_psi_2(t, out_param.d, out_param.radius, hyperbox))*volume;% the psi function is the transformation
        end
        hyperbox = bsxfun(@plus, zeros(2, out_param.d), [-inf; inf]);% the hyperbox must be the domain of the transformation, which is a this unit box
        out_param.measure = 'normal';% then a normal distribution can be used
    end
end

%------------------------------------------------------------------------------
% Minimum gathering of points
l_star = out_param.mmin - r_lag; % Minimum gathering of points for the sums of DFT
omg_circ = @(m) 2.^(-m);
omg_hat = @(m) out_param.fudge(m)/((1+out_param.fudge(r_lag))*omg_circ(r_lag));



FuncCount=cell2mat(FuncCount);
CVCount=cell2mat(CVCount);


% intialize CV param, redefine target function
mu=0;beta=0;
if cv.J % if using control variates(f is structure), redefine f
    mu = f.cv; f = f.func;
end

if strcmp(out_param.measure,'normal')
    f=@(x) f(gail.stdnorminv(x));
elseif strcmp(out_param.measure,'uniform')
    Cnorm = prod(hyperbox(2,:)-hyperbox(1,:));
    f=@(x) Cnorm*f(bsxfun(@plus,hyperbox(1,:),bsxfun(@times,(hyperbox(2,:)-hyperbox(1,:)),x))); % a + (b-a)x = u
end

if strcmp(out_param.transform,'Baker')
    f=@(x) f(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(out_param.transform,'C0')
    f=@(x) f(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(out_param.transform,'C1')
    f=@(x) f(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(out_param.transform,'C1sin')
    f=@(x) f(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
end

%% Main algorithm - Preallocation
Stilde=zeros(out_param.mmax-out_param.mmin+1,1); %initialize sum of DFT terms
CStilde_low = -inf(1,out_param.mmax-l_star+1); %initialize various sums of DFT terms for necessary conditions
CStilde_up = inf(1,out_param.mmax-l_star+1); %initialize various sums of DFT terms for necessary conditions
errest=zeros(out_param.mmax-out_param.mmin+1,1); %initialize error estimates
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1); %initialize approximations to integral
exit_len = 2;
out_param.exit=false(1,exit_len); %we start the algorithm with all warning flags down

%% Initial points and FFT
out_param.n=2^out_param.mmin; %total number of points to start with
n0=out_param.n; %initial number of points
xpts=mod(bsxfun(@plus, gail.lattice_gen(1,n0,out_param.d), out_param.shift),1); %grab Lattice points

y=f(xpts); %evaluate integrand
yval=y;


% evaluate integrand
if cv.J==0 % no control variates
    y = f(xpts); yval = y;
else  % using control variates
    ycv = f(xpts);
    y = ycv(:,1:FuncCount); yval = y;
    yg = ycv(:,FuncCount+1:end); %yvalg = yg;
end

%% Compute initial FFT
for l=0:out_param.mmin-1
    nl=2^l;
    nmminlm1=2^(out_param.mmin-l-1);
    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
    coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
    coefv=repmat(coef,nmminlm1,1);
    
    evenval=y(ptind);
    oddval=y(~ptind);
    y(ptind)=(evenval+coefv.*oddval)/2;
    y(~ptind)=(evenval-coefv.*oddval)/2;
    
    if cv.J
        evenval=yg(ptind, (1:cv.J));
        oddval=yg(~ptind, (1:cv.J));
        yg(ptind, (1:cv.J))=(evenval+coefv.*oddval)/2;
        yg(~ptind, (1:cv.J))=(evenval-coefv.*oddval)/2;
        
    end
    % y now contains the FFT coefficients
    
end

%% Create kappanumap implicitly from the data
kappanumap=(1:out_param.n)'; %initialize map

for l=out_param.mmin-1:-1:1
    nl=2^l;
    oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
    newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
    flip=find(newone>oldone); %which in the pair are the larger ones
    if ~isempty(flip)
        flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
        flipall=flipall(:);
        temp=kappanumap(nl+1+flipall); %then flip
        kappanumap(nl+1+flipall)=kappanumap(1+flipall); %them
        kappanumap(1+flipall)=temp; %around
    end
end

%% If using control variates, find optimal beta
if cv.J
    
    X = yg(kappanumap(2^(out_param.mmin-r_lag-1)+1:end), (1:end));
    Y =  y(kappanumap(2^(out_param.mmin-r_lag-1)+1:end), (1:end));
    
    Val = [];
    Val = [Y X];
    meanVal=[mean(Val)];
    
    A=bsxfun(@minus, Val, meanVal);
    
    C=[ones(FuncCount,1); zeros(CVCount,1)];
    [U,S,V]=svd([A; C'],0);
    Sdiag = diag(S);
    U2=U(end,:);
    H=U2'/(U2*U2');
    beta=V*(H./Sdiag);
    
    display(beta);
    
    meanX=meanVal(:,FuncCount+1:end);
    meanX=[zeros(FuncCount,1); meanX'];
    
    Ytemp=[];
    Ytemp=[y yg];
    
    yval = ycv(:,1:FuncCount)*beta(1:FuncCount,:) + ycv(:,FuncCount+1:end)*beta(FuncCount+1:end,:);
    y = (y(:,1:end)*beta(1:FuncCount,:)) + (yg(:,1:end)*beta(FuncCount+1:end,:));
    
    %% rebuild kappa map
    kappanumap=(1:out_param.n); %initialize map
    
    for l=out_param.mmin-1:-1:1
        
        nl=2^l;
        oldone=abs(y(kappanumap(2:nl)));
        newone=abs(y(kappanumap(nl+2:2*nl)));
        
        flip=find(newone>oldone);
        if ~isempty(flip)
            flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
            flipall=flipall(:);
            temp=kappanumap(nl+1+flipall); %then flip
            kappanumap(nl+1+flipall)=kappanumap(1+flipall); %them
            kappanumap(1+flipall)=temp; %around
        end
    end
end

%% Compute Stilde (1)
nllstart=int64(2^(out_param.mmin-r_lag-1));
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
out_param.bound_err=out_param.fudge(out_param.mmin)*Stilde(1);
errest(1)=out_param.bound_err;

% Necessary conditions
for l = l_star:out_param.mmin % Storing the information for the necessary conditions
    C_low = 1/(1+omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    C_up = 1/(1-omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    CStilde_low(l-l_star+1) = max(CStilde_low(l-l_star+1),C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
    if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
        CStilde_up(l-l_star+1) = min(CStilde_up(l-l_star+1),C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
    end
end
if any(CStilde_low(:) > CStilde_up(:))
    out_param.exit(2) = true;
end

%% Approximate integral (1)
display('approximate integral (1)');
if cv.J
    q= mean(yval) - mu*beta(FuncCount+1:end,:);
else
    q=mean(yval);
end

% Check the end of the algorithm
q = q - errest(1)*(max(out_param.abstol, out_param.reltol*abs(q + errest(1)))...
    - max(out_param.abstol, out_param.reltol*abs(q - errest(1))))/...
    (max(out_param.abstol, out_param.reltol*abs(q + errest(1)))...
    + max(out_param.abstol, out_param.reltol*abs(q - errest(1)))); % Optimal estimator


display(q);
%q=q(1);
appxinteg(1)=q;

is_done = false;
if 4*errest(1)^2/(max(out_param.abstol, out_param.reltol*abs(q + errest(1)))...
        + max(out_param.abstol, out_param.reltol*abs(q - errest(1))))^2 <= 1
    out_param.time=toc(t_start);
    is_done = true;
elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
    out_param.exit(1) = true;
    is_done = true;
end

%% Loop over m
for m=out_param.mmin+1:out_param.mmax
    if is_done,
        break;
    end
    
    out_param.n=2^m;
    mnext=m-1;
    nnext=2^mnext;
    xnext=mod(bsxfun(@plus, gail.lattice_gen(nnext+1,2*nnext,out_param.d), out_param.shift),1);
    n0=n0+nnext;
    
    % check for using control variates or not
    if cv.J == 0
        ynext = f(xnext); yval=[yval; ynext];
    else
        
        ycvnext = f(xnext);
        ynext = ycvnext(:,1:FuncCount)*beta(1:FuncCount,:) + ycvnext(:,FuncCount+1:end)*beta(FuncCount+1:end,:);
        yval=[yval; ynext];
        
        
        
    end
    
    %% Compute initial FFT on next points
    % Not updating beta
    if out_param.betaUpdate==0
        for l=0:mnext-1
            nl=2^l;
            nmminlm1=2^(mnext-l-1);
            ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
            coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
            coefv=repmat(coef,nmminlm1,1);
            evenval=ynext(ptind);
            oddval=ynext(~ptind);
            ynext(ptind)=(evenval+coefv.*oddval)/2;
            ynext(~ptind)=(evenval-coefv.*oddval)/2;
        end
        
        %Compute FFT on all points
        y=[y;ynext];
        nl=2^mnext;
        ptind=[true(nl,1); false(nl,1)];
        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
        coefv=repmat(coef,nmminlm1,1);
        evenval=y(ptind);
        oddval=y(~ptind);
        y(ptind)=(evenval+coefv.*oddval)/2;
        y(~ptind)=(evenval-coefv.*oddval)/2;
        
        % Update kappanumap
        kappanumap=[kappanumap; 2^(m-1)+kappanumap]; %initialize map
        for l=m-1:-1:m-r_lag
            nl=2^l;
            oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
            newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
            flip=find(newone>oldone);
            if ~isempty(flip)
                flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
                flipall=flipall(:);
                temp=kappanumap(nl+1+flipall);
                kappanumap(nl+1+flipall)=kappanumap(1+flipall);
                kappanumap(1+flipall)=temp;
            end
        end
        
        % Updating beta
    else
        ycv = [ycv;ycvnext];y = ycv(:,1);yg = ycv(:,2:end);
        
        % compute FFT
        y=[y;ynext];
        nl=2^mnext;
        ptind=[true(nl,1); false(nl,1)];
        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
        coefv=repmat(coef,nmminlm1,1);
        evenval=y(ptind);
        oddval=y(~ptind);
        y(ptind)=(evenval+coefv.*oddval)/2;
        y(~ptind)=(evenval-coefv.*oddval)/2;
        
        for l=0:m-1
            nl=2^l;
            nmminlm1=2^(m-l-1);
            ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
            evenval=y(ptind);
            oddval=y(~ptind);
            coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
            coefv=repmat(coef,nmminlm1,1);
            y(ptind)=(evenval+coefv.*oddval)/2;
            y(~ptind)=(evenval-coefv.*oddval)/2;
            evenval=yg(ptind, (1:cv.J));
            oddval=yg(~ptind, (1:cv.J));
            yg(ptind, (1:cv.J))=(evenval+coefv.*oddval)/2;
            yg(~ptind, (1:cv.J))=(evenval-coefv.*oddval)/2;
        end
        
        X = yg(kappanumap(2^(m-r_lag-1)+1:end), (1:cv.J));
        Y = y(kappanumap(2^(m-r_lag-1)+1:end));
        beta = real(X \ Y);
        out_param.beta = [out_param.beta;beta];
        yval = ycv(:,1) - ycv(:,2:end)*beta;
        y = y-yg*beta;
        
        %% update kappamap
        kappanumap=[kappanumap; 2^(m-1)+kappanumap]; %initialize map
        for l=m-1:-1:m-r_lag
            nl=2^l;
            oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
            newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
            flip=find(newone>oldone);
            if ~isempty(flip)
                flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
                flipall=flipall(:);
                temp=kappanumap(nl+1+flipall);
                kappanumap(nl+1+flipall)=kappanumap(1+flipall);
                kappanumap(1+flipall)=temp;
            end
        end
        
    end
    
    %% Compute Stilde (2)
    nllstart=int64(2^(m-r_lag-1));
    meff=m-out_param.mmin+1;
    Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
    out_param.bound_err=out_param.fudge(m)*Stilde(meff);
    errest(meff)=out_param.bound_err;
    
    % Necessary conditions
    for l = l_star:m % Storing the information for the necessary conditions
        C_low = 1/(1+omg_hat(m-l)*omg_circ(m-l));
        C_up = 1/(1-omg_hat(m-l)*omg_circ(m-l));
        CStilde_low(l-l_star+1) = max(CStilde_low(l-l_star+1),C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
        if (omg_hat(m-l)*omg_circ(m-l) < 1)
            CStilde_up(l-l_star+1) = min(CStilde_up(l-l_star+1),C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
        end
    end
    
    if any(CStilde_low(:) > CStilde_up(:))
        out_param.exit(2) = true;
    end
    
    %% Approximate integral (2)
    if cv.J
        q= mean(yval) - mu*beta(FuncCount+1:end,:);
    else
        q=mean(yval);
    end
    % Check the end of the algorithm
    q = q - errest(meff)*(max(out_param.abstol, out_param.reltol*abs(q + errest(meff)))...
        - max(out_param.abstol, out_param.reltol*abs(q - errest(meff))))/...
        (max(out_param.abstol, out_param.reltol*abs(q + errest(meff)))...
        + max(out_param.abstol, out_param.reltol*abs(q - errest(meff)))); % Optimal estimator
    
    appxinteg(meff)=q;
    
    if 4*errest(meff)^2/(max(out_param.abstol, out_param.reltol*abs(q + errest(meff)))...
            + max(out_param.abstol, out_param.reltol*abs(q - errest(meff))))^2 <= 1
        out_param.time=toc(t_start);
        is_done = true;
    elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
        out_param.exit(1) = true;
    end
end

% Decode the exit structure
exit_str=2.^(0:exit_len-1).*out_param.exit;
exit_str(out_param.exit==0)=[];
if numel(exit_str)==0;
    out_param.exitflag=0;
else
    out_param.exitflag=exit_str;
end

out_param = rmfield(out_param,'exit');
out_param.time=toc(t_start);

end


%% Parsing for the input of cubLattice_g
function [f,hyperbox, out_param,cv, FuncCount, CVCount] = cubLattice_g_param(r_lag,varargin)

FuncCount=varargin(6);
CVCount=varargin(7);

% Default parameter values
default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
default.measure  = 'uniform';
default.transf = 1;% default transformation (box-to-ball or box-to-sphere)
default.radius = 1;% radius of the ball or sphere
default.abstol  = 1e-4;
default.reltol  = 1e-2;
default.shift  = rand;
default.mmin  = 10;
default.mmax  = 20;
default.fudge = @(m) 10*2.^-(m);
default.transform = 'Baker';
default.betaUpdate=0;

% two data structures for function: function || structure(using CV)
validf = @(x) gail.isfcn(x) || isstruct(x);

if numel(varargin)<2
    help cubLattice_g
    warning('GAIL:cubLattice_g:fdnotgiven',...
        'At least, function f and hyperbox need to be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    hyperbox = default.hyperbox;
else
    f = varargin{1};
    if ~validf(f)
        warning('GAIL:cubLattice_g:fnotfcn',...
            'The given input f should be a function, or a structure in required form.' )
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
    elseif numel(varargin) == 2
        out_param.f=f;
        hyperbox = varargin{2};
        if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(size(hyperbox,2)<601)
            warning('GAIL:cubLattice_g:hyperbox_error1',...
                'The hyperbox must be a real matrix of size 2xd where d can not be greater than 600. Example for f(x)=x^2:')
            f = @(x) x.^2;
            out_param.f=f;
            hyperbox = default.hyperbox;
        end
    else
        out_param.f = f;
        hyperbox = varargin{2}; % hyperbox validation will be done above
        default.shift = rand(1, size(hyperbox, 2)); % The shift needs to take the dimension of the problem
    end
end

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin(3:end);
    for j=1:numel(varargin)-2
        validvarargin=validvarargin && (isnumeric(in3{j}) ...
            || ischar(in3{j}) || isstruct(in3{j}) || gail.isfcn(in3{j}));
    end
    if ~validvarargin
        warning('GAIL:cubLattice_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in3=varargin{3};
end

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end

if ~validvarargin
    out_param.measure = default.measure;
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.shift = default.shift;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;
    out_param.fudge = default.fudge;
    out_param.transform = default.transform;
    out_param.betaUpdate=default.betaUpdate;
else
    p = inputParser;
    addRequired(p,'f', validf);
    addRequired(p,'hyperbox',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal','uniform ball',...
            'uniform ball_box','uniform ball_normal','uniform sphere',...
            'uniform sphere_box','uniform sphere_normal'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'shift',default.shift,@isnumeric);
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
        addOptional(p,'transform',default.transform,...
            @(x) any(validatestring(x, {'id','Baker','C0','C1','C1sin'})));
        addOptional(p,'betaUpdate',default.betaUpdate,@isnumeric);
        
    else
        if isstruct(in3) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal','uniform ball',...
            'uniform ball_box','uniform ball_normal','uniform sphere',...
            'uniform sphere_box','uniform sphere_normal'})));
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'shift',default.shift,@isnumeric);
        f_addParamVal(p,'mmin',default.mmin,@isnumeric);
        f_addParamVal(p,'mmax',default.mmax,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@gail.isfcn);
        f_addParamVal(p,'transform',default.transform,...
            @(x) any(validatestring(x, {'id','Baker','C0','C1','C1sin'})));
        f_addParamVal(p,'betaUpdate',default.betaUpdate,@isnumeric);
        
    end
    parse(p,f,hyperbox,varargin{3:end});
    out_param = p.Results;
end

% get the number of control variates
if ~isstruct(f) %  not using control variates
    cv.J = 0;
else % using control variates, checking mu
    if isnumeric(f.cv)
        cv.J = size(f.cv,2);
    else
        warning('GAIL:cubLattice_g:controlvariates_error1',...
            'f.cv should be numerical values');
    end
end

% Force measure to be one of the allowed ones
if ~(strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') || ...
        strcmp(out_param.measure,'uniform ball') || ...
        strcmp(out_param.measure,'uniform ball_box') || ...
        strcmp(out_param.measure,'uniform ball_normal') || ...
        strcmp(out_param.measure,'uniform sphere') || ...
        strcmp(out_param.measure,'uniform sphere_box') || ...
        strcmp(out_param.measure,'uniform sphere_normal'))
    warning('GAIL:cubLattice_g:notmeasure',['Given measure is not allowed.' ...
        ' Using default measure ' num2str(default.measure)])
    out_param.measure = default.measure;
end

% simplifying out_param.measure and storing which transformation should be
% applied
if strcmp(out_param.measure,'uniform ball') || strcmp(out_param.measure,'uniform sphere')
    out_param.transf = default.transf;
elseif strcmp(out_param.measure,'uniform ball_box')
    out_param.transf = 1; % one means from a box
    out_param.measure = 'uniform ball';
elseif strcmp(out_param.measure,'uniform ball_normal')
    out_param.transf = 2; % two means from normal
    out_param.measure = 'uniform ball';
elseif strcmp(out_param.measure,'uniform sphere_box')
    out_param.transf = 1;
    out_param.measure = 'uniform sphere';
elseif strcmp(out_param.measure,'uniform sphere_normal')
    out_param.transf = 2;
    out_param.measure = 'uniform sphere';
end

%validating hyperbox
if strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') % 'uniform box' or 'normal box'
    out_param.d = size(hyperbox,2);
    if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(out_param.d<1111)
        warning('GAIL:cubLattice_g:hyperbox_error2',...
            'When measure is ''uniform'' or ''normal'', the hyperbox must be a real matrix of size 2 x d where d can not be greater than 1111. Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end
else % 'uniform ball' or 'uniform sphere'
    out_param.d = size(hyperbox,2);
    if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==1) || ~(out_param.d<1112) % size(hyperbox,2) is actually equal to d+1 (the extra value is the radius)
        warning('GAIL:cubLattice_g:hyperbox_error3',...
            'When measure is ''uniform ball'' or ''uniform sphere'', the hyperbox must be a real matrix of size 1 x (d+1) where d can not be greater than 1111.  Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end
    
    out_param.radius = hyperbox(out_param.d);
    hyperbox = hyperbox(:,1:out_param.d-1); % removing the last value is the radius, which is the radius
    out_param.d = out_param.d - 1; % storing the rigth dimension of the ball or sphere
    out_param.shift = out_param.shift(1:end-1);
    
    if strcmp(out_param.measure,'uniform ball') && out_param.d <= 0
        warning('GAIL:cubLattice_g:dimensionequalszero',...
            'When measure is ''uniform ball'', the number of dimentions must be a positive values.  Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end
    
    if strcmp(out_param.measure,'uniform sphere') && out_param.d <= 1
        warning('GAIL:cubLattice_g:dimensionsmallerthan2',...
            'When measure is ''uniform sphere'', the number of dimentions must be at least 2.  Example for f(x)=x^2 over [0,1]:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
        out_param.d = size(hyperbox,2);
        out_param.measure = default.measure;
    end
    
    if ~isfinite(out_param.radius) || out_param.radius <= 0.0
        warning('GAIL:cubLattice_g:infiniteradius',...
            'When measure is ''uniform ball'' or ''uniform sphere'', the radius must a finite positive real number. Default value for the radius will be used:')
        out_param.radius = default.radius;
    end
    
    if strcmp(out_param.measure,'uniform sphere') && out_param.transf == 1 % box-to-sphere transformation
        % setting out_param.d to be the dimension of the box over which the
        % integral will actually be computed
        out_param.d = out_param.d - 1;
        out_param.shift = out_param.shift(1:end-1);
    end
end

% Force absolute tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:cubLattice_g:abstolnonpos',['Absolute tolerance cannot be negative.' ...
        ' Using default absolute tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol < 0) || (out_param.reltol > 1)
    warning('GAIL:cubLattice_g:reltolnonunit',['Relative tolerance should be chosen in [0,1].' ...
        ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Checks if shift is a vector in [0,1)^d
if ~(all(out_param.shift < 1) && all(out_param.shift >= 0) && size(out_param.shift, 2) == out_param.d)
    warning('GAIL:cubLattice_g:shift',['The shift should be a vector ' ...
        'of size 1 x d in [0,1)^d.' ...
        ' Using default shift ' num2str(default.shift)])
    out_param.shift = default.shift;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('GAIL:cubLattice_g:lowmmin',['The minimum starting exponent ' ...
        'should be an integer greater than 0 and smaller or equal than the maxium.' ...
        ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force mmin to be integer greater than r_lag (so that l_star=mmin-r_lag>=0)
if out_param.mmin < r_lag
    warning('GAIL:cubLattice_g:lowmminrlag',['The minimum starting exponent ' ...
        'should be at least ' num2str(r_lag) '.' ...
        ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 20
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin)
    warning('GAIL:cubLattice_g:wrongmmax',['The maximum exponent for the budget should be an integer bigger than mmin and smaller than the allowed by gail.lattice_gen.' ...
        ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('GAIL:cubLattice_g:fudgenonpos',['The fudge factor should be a positive function.' ...
        ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force transform to only be id, Baker, C0, C1 or C1sin
if ~(strcmp(out_param.transform,'id') || strcmp(out_param.transform,'Baker') || strcmp(out_param.transform,'C0') || strcmp(out_param.transform,'C1') || strcmp(out_param.transform,'C1sin') )
    warning('GAIL:cubLattice_g:notmeasure',['The periodizing transformations can only be id, Baker, C0, C1 or C1sin.' ...
        ' Using default error tolerance ' num2str(default.transform)])
    out_param.transform = default.transform;
end

% Checking on pure absolute/relative error
if (out_param.abstol==0) && (out_param.reltol==0)
    warning('GAIL:cubLattice_g:tolzeros',['Absolute and relative error tolerances can not be simultaniusly 0.' ...
        ' Using default absolute tolerance ' num2str(default.abstol) ' and relative tolerance ' num2str(default.reltol)])
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
end

% Checking on the hyperbox given the measure
if (strcmp(out_param.measure,'uniform')) && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubLattice_g:hyperboxnotfinite',['If uniform measure, hyperbox must be of finite volume.' ...
        ' Using default hyperbox:'])
    disp([zeros(1,out_param.d);ones(1,out_param.d)])
    hyperbox = [zeros(1,out_param.d);ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(any(isfinite(hyperbox)))>0)
    warning('GAIL:cubLattice_g:hyperboxfinite',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
        ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(hyperbox(1,:)==hyperbox(2,:)) || any(hyperbox(1,:)>hyperbox(2,:)))
    warning('GAIL:cubLattice_g:hyperboxnormalwrong',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
        ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'uniform ball') || strcmp(out_param.measure,'uniform sphere'))...
        && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubLattice_g:infinitecoordinateforthecenter',['If uniform ball or sphere measure, all the coordinates of the center must be finite.' ...
        ' Using the origin as the center:'])
    % out_param.d should not be used here because this variable stores the
    % dimension of the box over which the integral will actually be
    % computed, whih may be different from the dimesion of the sphere
    hyperbox = zeros(1,size(hyperbox,2));
end
end
