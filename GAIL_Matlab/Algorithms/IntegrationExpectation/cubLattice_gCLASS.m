function [meanf, mean_out] = cubLattice_gCLASS(varargin)
% Example 1: 
% >> w.f= @(x) prod(x,2); w.absTol=1e-5;
% >> w.relTol=0; w.domain = [zeros(1,2);ones(1,2)];
% >> w.transform = 'C1sin'; 
% >> [q, out_param] = cubLattice_gCLASS(w);
% >> exactsol = 1/4;
% >> check = abs(exactsol-q) < 1e-5
% check = 1
%
% Example 2: 
% >> w.f= @(x) exp(-x(:,1).^2-x(:,2).^2); w.absTol=1e-3;
% >> w.relTol=1e-2; w.domain = [-ones(1,2);2*ones(1,2)];
% >> w.transform = 'C1';
% >> [q, out_param] = cubLattice_gCLASS(w);
% >> exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
% >> check = abs(exactsol-q) < 1e-3;
% check = 1 
%
% Example 3: 
% >> w.f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); w.domain = [-inf(1,1);inf(1,1)];
% >> w.measure='normal'; w.absTol=1e-4; w.relTol=1e-2;
% w.transform=('C1sin');
% >> [q, out_param] = cubLattice_gCLASS(w);
% >> price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
% >> check = abs(price-q) < 1e-4;
% check = 1
%
% Example 4: 
% >> w.f = @(x) 8*prod(x,2); w.domain = [zeros(1,5);ones(1,5)];
% >> w.measure='uniform'; 
% >> w.absTol=1e-5; w.relTol=0;
% >> [q, out_param] = cubLattice_gCLASS(w);
% >> exactsol = 1/4;
% >> check = abs(exactsol-q) < 1e-5;
% check = 1
%
% Example 5: 
% >> w.f= @(x) 3./(5-4*(cos(2*pi*x)));
% >> w.absTol=1e-5; w.relTol=0;
% >> w.domain = [0;1]; w.transform = 'id';
% >> [q, out_param] = cubLattice_gCLASS(w);
% >> exactsol = 1;
% >> check = abs(exactsol-q) < 1e-5;
% check = 1
%
% Example 6: 
% >> w.f= @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2]; w.domain = [zeros(1,3);2*ones(1,3)];
% >> w.trueMuCV=[8,32/3]; w.absTol=1e-6;
% >> w.relTol=0; w.measure = 'uniform';
% >> [q, out_param] = cubLattice_gCLASS(w);
% >> exactsol = 128/3;
% >> check = abs(exactsol-q) < 1e-6;
% check = 1

t_start = tic;
%% Initial important cone factors and Check-initialize parameters
r_lag = 4; %distance between coefficients summed and those computed

mean_inp = gail.cubLatticeParam(varargin{:}); %parse the input and check it for errors
mean_inp.fun.nMax = min(mean_inp.fun.nMax,2^24);
mean_out = gail.cubLatticeOut(mean_inp); %create the output class

%------------------------------------------------------------------------------
% Minimum gathering of points
l_star = mean_out.mmin - r_lag; % Minimum gathering of points for the sums of DFT
omg_circ = @(m) 2.^(-m);
omg_hat = @(m) mean_out.CM.inflateFun(m)/((1+mean_out.CM.inflateFun(r_lag))*omg_circ(r_lag));

% intialize CV param, redefine target function
mu=0;beta=0;

if mean_out.CM.nCV  % if using control variates(f is structure), redefine f
    mu = mean_out.CM.trueMuCV;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main algorithm - Preallocation
Stilde=zeros(mean_out.mmax-mean_out.mmin+1,1); %initialize sum of DFT terms
CStilde_low = -inf(1,mean_out.mmax -l_star+1); %initialize various sums of DFT terms for necessary conditions
CStilde_up = inf(1,mean_out.mmax -l_star+1); %initialize various sums of DFT terms for necessary conditions
errest=zeros(mean_out.mmax -mean_out.mmin+1,1); %initialize error estimates
appxinteg=zeros(mean_out.mmax -mean_out.mmin+1,1); %initialize approximations to integral
exit_len = 2;
exit=false(1,exit_len); %we start the algorithm with all warning flags down


%% Initial points and FFT
mean_out.nSample=2^mean_out.mmin; %total number of points to start with
n0=mean_out.nSample; %initial number of points
xpts=mod(bsxfun(@plus, gail.lattice_gen(1,n0,mean_out.fun.d), mean_out.shiftVal),1); %grab Lattice points

y=mean_out.fff(xpts); %evaluate integrand
yval=y;

% evaluate integrand
ycv = mean_out.fff(xpts);
y = ycv(:,1:mean_out.nf); yval = y;
yg = ycv(:,mean_out.nf+1:end); %yvalg = yg;


%% Compute initial FFT
for l=0:mean_out.mmin-1
    nl=2^l;
    nmminlm1=2^(mean_out.mmin-l-1);
    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
    coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
    coefv=repmat(coef,nmminlm1,1);
    
    evenval=y(ptind);
    oddval=y(~ptind);
    y(ptind)=(evenval+coefv.*oddval)/2;
    y(~ptind)=(evenval-coefv.*oddval)/2;
    
    if mean_out.CM.nCV 
        evenval=yg(ptind, (1:mean_out.CM.nCV ));
        oddval=yg(~ptind, (1:mean_out.CM.nCV ));
        yg(ptind, (1:mean_out.CM.nCV ))=(evenval+coefv.*oddval)/2;
        yg(~ptind, (1:mean_out.CM.nCV ))=(evenval-coefv.*oddval)/2;
        
    end
    % y now contains the FFT coefficients
end

%% Create kappanumap implicitly from the data
kappanumap=(1:mean_out.nSample)'; %initialize map

for l=mean_out.mmin-1:-1:1
    nl=2^l;
    oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
    newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa,
    flip=find(newone>oldone); %which in the pair are the larger ones
    if ~isempty(flip)
        flipall=bsxfun(@plus,flip,0:2^(l+1):2^mean_out.mmin-1);
        flipall=flipall(:);
        temp=kappanumap(nl+1+flipall); %then flip
        kappanumap(nl+1+flipall)=kappanumap(1+flipall); %them
        kappanumap(1+flipall)=temp; %around
    end
end

%% Finding optimal beta
% Pre-determine the size of the beta coefficients 
if mean_out.CM.nCV 
    C=[ones(mean_out.nf,1); zeros(mean_out.CM.nCV,1)];
else 
    C=[ones(mean_out.nf,1)];
end 
    
%% alogirhtm to find beta 
X = yg(kappanumap(2^(mean_out.mmin-r_lag-1)+1:end), (1:end));
Y =  y(kappanumap(2^(mean_out.mmin-r_lag-1)+1:end), (1:end));

Val = [];
Val = [Y X];
meanVal=[mean(Val)];

A=bsxfun(@minus, Val, meanVal);

[U,S,V]=svd([A; C'],0);
Sdiag = diag(S);
U2=U(end,:);
H=U2'/(U2*U2');
beta=V*(H./Sdiag);
beta=real(beta);

meanX=meanVal(:,mean_out.nf+1:end);
meanX=[zeros(mean_out.nf,1); meanX'];

Ytemp=[];
Ytemp=[y yg];

yval = ycv(:,1:mean_out.nf)*beta(1:mean_out.nf,:) + ycv(:,mean_out.nf+1:end)*beta(mean_out.nf+1:end,:);
y = (y(:,1:end)*beta(1:mean_out.nf,:)) + (yg(:,1:end)*beta(mean_out.nf+1:end,:));

%% rebuild kappa map
kappanumap=(1:mean_out.nSample); %initialize map

for l=mean_out.mmin-1:-1:1
    
    nl=2^l;
    oldone=abs(y(kappanumap(2:nl)));
    newone=abs(y(kappanumap(nl+2:2*nl)));
    
    flip=find(newone>oldone);
    if ~isempty(flip)
        flipall=bsxfun(@plus,flip,0:2^(l+1):2^mean_out.mmin-1);
        flipall=flipall(:);
        temp=kappanumap(nl+1+flipall); %then flip
        kappanumap(nl+1+flipall)=kappanumap(1+flipall); %them
        kappanumap(1+flipall)=temp; %around
    end
end


%% Compute Stilde (1)
nllstart=int64(2^(mean_out.mmin-r_lag-1));
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
mean_out.errBd=mean_out.CM.inflateFun(mean_out.mmin)*Stilde(1);

errest(1)=mean_out.errBd;

% Necessary conditions
for l = l_star:mean_out.mmin % Storing the information for the necessary conditions
    C_low = 1/(1+omg_hat(mean_out.mmin-l)*omg_circ(mean_out.mmin-l));
    C_up = 1/(1-omg_hat(mean_out.mmin-l)*omg_circ(mean_out.mmin-l));
    CStilde_low(l-l_star+1) = max(CStilde_low(l-l_star+1),C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
    if (omg_hat(mean_out.mmin-l)*omg_circ(mean_out.mmin-l) < 1)
        CStilde_up(l-l_star+1) = min(CStilde_up(l-l_star+1),C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l)))));
    end
end
if any(CStilde_low(:) > CStilde_up(:))
    exit(2) = true;
end

%% Approximate integral (1)
if mean_out.CM.nCV 
    q = mean(yval) - mu*beta(mean_out.nf+1:end,:);
else
    q =mean(yval);
end

% Check the end of the algorithm
q = q - errest(1)*(max(mean_out.err.absTol, mean_out.err.relTol*abs(q + errest(1)))...
    - max(mean_out.err.absTol, mean_out.err.relTol*abs(q - errest(1))))/...
    (max(mean_out.err.absTol, mean_out.err.relTol*abs(q + errest(1)))...
    + max(mean_out.err.absTol, mean_out.err.relTol*abs(q - errest(1)))); % Optimal estimator

q=q(1);
appxinteg(1)=q;

is_done = false;
if 4*errest(1)^2/(max(mean_out.err.absTol, mean_out.err.relTol*abs(q + errest(1)))...
        + max(mean_out.err.absTol, mean_out.err.relTol*abs(q - errest(1))))^2 <= 1
    mean_out.time=toc(t_start);
    is_done = true;
elseif mean_out.mmin == mean_out.mmax  % We are on our max budget and did not meet the error condition => overbudget
    exit(1) = true;
    is_done = true;
end

%% Loop over m
for m=mean_out.mmin+1:mean_out.mmax 
    if is_done
        break;
    end
    
    mean_out.nSample=2^m;
    mnext=m-1;
    nnext=2^mnext;
    xnext=mod(bsxfun(@plus, gail.lattice_gen(nnext+1,2*nnext,mean_out.fun.d), mean_out.shiftVal),1);
    n0=n0+nnext;
    
    % check for using control variates or not
    if mean_out.CM.nCV  == 0
        % ynext = f(xnext); yval=[yval; ynext];
        ycvnext = mean_out.fff(xnext);
        ynext = ycvnext(:,1:mean_out.nf)*beta(1:mean_out.nf,:) + ycvnext(:,mean_out.nf+1:end)*beta(mean_out.nf+1:end,:);
        yval=[yval; ynext];
    else
        ycvnext = mean_out.fff(xnext);
        ynext = ycvnext(:,1:mean_out.nf)*beta(1:mean_out.nf,:) + ycvnext(:,mean_out.nf+1:end)*beta(mean_out.nf+1:end,:);
        yval=[yval; ynext];
    end
    
    %% Compute initial FFT on next points
    % Not updating beta
    if mean_out.betaUpdate==0
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
            evenval=yg(ptind, (1:mean_out.CM.nCV ));
            oddval=yg(~ptind, (1:mean_out.CM.nCV ));
            yg(ptind, (1:mean_out.CM.nCV ))=(evenval+coefv.*oddval)/2;
            yg(~ptind, (1:mean_out.CM.nCV ))=(evenval-coefv.*oddval)/2;
        end
        
        X = yg(kappanumap(2^(m-r_lag-1)+1:end), (1:mean_out.CM.nCV ));
        Y = y(kappanumap(2^(m-r_lag-1)+1:end));
        beta = real(X \ Y);
        mean_out.sol = [mean_out.sol;beta];
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
    meff=m-mean_out.mmin+1;
    Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
    mean_out.errBd=mean_out.CM.inflateFun(m)*Stilde(meff);
    errest(meff)=mean_out.errBd;
    
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
        exit(2) = true;
    end
    
    %% Approximate integral (2)
    if mean_out.CM.nCV 
        q= mean(yval) - mu*beta(mean_out.nf+1:end,:);
    else
        q=mean(yval);
    end
    
    % Check the end of the algorithm
    q = q - errest(meff)*(max(mean_out.err.absTol, mean_out.err.relTol*abs(q + errest(meff)))...
        - max(mean_out.err.absTol, mean_out.err.relTol*abs(q - errest(meff))))/...
        (max(mean_out.err.absTol, mean_out.err.relTol*abs(q + errest(meff)))...
        + max(mean_out.err.absTol, mean_out.err.relTol*abs(q - errest(meff)))); % Optimal estimator
    
    appxinteg(meff)=q;
    
    if 4*errest(meff)^2/(max(mean_out.err.absTol, mean_out.err.relTol*abs(q + errest(meff)))...
            + max(mean_out.err.absTol, mean_out.err.relTol*abs(q - errest(meff))))^2 <= 1
        mean_out.time=toc(t_start);
        is_done = true;
    elseif m == mean_out.mmax  % We are on our max budget and did not meet the error condition => overbudget
        exit(1) = true;
    end
end

% Decode the exit structure
exit_str=2.^(0:exit_len-1).*exit;
exit_str(exit==0)=[];
if numel(exit_str)==0;
    exitflag=0;
else
    exitflag=exit_str;
end

mean_out.time=toc(t_start);
mean_out.sol=q;

meanf=mean_out.sol;

end


