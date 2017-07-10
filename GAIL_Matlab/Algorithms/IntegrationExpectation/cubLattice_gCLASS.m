function [meanf, mean_out] = cubLattice_gCLASS(varargin)

t_start = tic;
%% Initial important cone factors and Check-initialize parameters
r_lag = 4; %distance between coefficients summed and those computed

mean_inp = gail.cubLatticeParam(varargin{:}); %parse the input and check it for errors
mean_out = gail.cubLatticeOut(mean_inp); %create the output class

%------------------------------------------------------------------------------
% % TRANSFORMATION
% % changing the integrand and the mean_out.domain when measure is uniform ball or
% % sphere by applying the appropriate transformation
% if strcmpi(mean_out.measure,'uniform ball') || strcmpi(mean_out.measure,'uniform sphere')% using uniformly distributed samples on a ball or sphere
%     if strcmp(mean_out.measure,'uniform sphere') && mean_out.transf == 1 %box-to-sphere transformation
%         mean_out.d = mean_out.d + 1; % changing mean_out.d to the dimension of the sphere
%         mean_out.shiftVal = [mean_out.shiftVal rand];
%     end
%     
%     if strcmpi(mean_out.measure,'uniform ball')% using the formula of the volume of a ball
%         volume = ((2.0*pi^(mean_out.d/2.0))/(mean_out.d*gamma(mean_out.d/2.0)))*mean_out.radius^mean_out.d; %volume of a d-dimentional ball
%     else % using the formula of the volume of a sphere
%         volume = ((2.0*pi^(mean_out.d/2.0))/(gamma(mean_out.d/2.0)))*mean_out.radius^(mean_out.d - 1); %volume of a d-dimentional sphere
%     end
%     
%     if mean_out.transf == 1 % box-to-ball or box-to-sphere transformation should be used
%         if mean_out.d == 1 % It is not necessary to multiply the function f by the volume, since no transformation is being made
%             mean_out.domain = [mean_out.domain - mean_out.radius; mean_out.domain + mean_out.radius];% for one dimension, the ball is actually an interval
%             mean_out.measure = 'uniform';% then a uniform distribution on a box can be used
%         else
%             if strcmpi(mean_out.measure,'uniform ball') % box-to-ball transformation
%                 f = @(t) f(gail.domain_balls_spheres.ball_psi_1(t, mean_out.d, mean_out.radius, mean_out.domain))*volume;% the psi function is the transformation
%             else %  % box-to-sphere transformation
%                 f = @(t) f(gail.domain_balls_spheres.sphere_psi_1(t, mean_out.d, mean_out.radius, mean_out.domain))*volume;% the psi function is the transformation
%                 mean_out.d = mean_out.d - 1;% the box-to-sphere transformation takes points from a (d-1)-dimensional box to a d-dimensional sphere
%                 mean_out.shiftVal = mean_out.shiftVal(1:end-1);
%             end
%             mean_out.domain = [zeros(1, mean_out.d); ones(1, mean_out.d)];% the mean_out.domain must be the domain of the transformation, which is a unit box
%             mean_out.measure = 'uniform';% then a uniform distribution on a box can be used
%         end
%     else % normal-to-ball or normal-to-sphere transformation should be used
%         if strcmpi(mean_out.measure,'uniform ball') % normal-to-ball transformation
%             f = @(t) f(gail.domain_balls_spheres.ball_psi_2(t, mean_out.d, mean_out.radius, mean_out.domain))*volume;% the psi function is the transformation
%         else % normal-to-sphere transformation
%             f = @(t) f(gail.domain_balls_spheres.sphere_psi_2(t, mean_out.d, mean_out.radius, mean_out.domain))*volume;% the psi function is the transformation
%         end
%         mean_out.domain = bsxfun(@plus, zeros(2, mean_out.d), [-inf; inf]);% the mean_out.domain must be the domain of the transformation, which is a this unit box
%         mean_out.measure = 'normal';% then a normal distribution can be used
%     end
% end

%------------------------------------------------------------------------------
% Minimum gathering of points
l_star = mean_out.mmin - r_lag; % Minimum gathering of points for the sums of DFT
omg_circ = @(m) 2.^(-m);
omg_hat = @(m) mean_out.inflate(m)/((1+mean_out.inflate(r_lag))*omg_circ(r_lag));

% intialize CV param, redefine target function
mu=0;beta=0;

if mean_out.nCV  % if using control variates(f is structure), redefine f
    mu = mean_out.trueMuCV;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALREADY INCLUDED 
% if strcmp(mean_out.measure,'normal')
%     f=@(x) f(gail.stdnorminv(x));
% elseif strcmp(mean_out.measure,'uniform')
%     Cnorm = prod(mean_out.domain(2,:)-mean_out.domain(1,:));
%     f=@(x) Cnorm*f(bsxfun(@plus,mean_out.domain(1,:),bsxfun(@times,(mean_out.domain(2,:)-mean_out.domain(1,:)),x))); % a + (b-a)x = u
% end

% if strcmp(mean_out.periodTransform,'Baker')
%     f=@(x) f(1-2*abs(x-1/2)); % Baker's transform
% elseif strcmp(mean_out.periodTransform,'C0')
%     f=@(x) f(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
% elseif strcmp(mean_out.periodTransform,'C1')
%     f=@(x) f(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
% elseif strcmp(mean_out.periodTransform,'C1sin')
%     f=@(x) f(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
% end

%% Main algorithm - Preallocation
Stilde=zeros(mean_out.mmax -mean_out.mmin+1,1); %initialize sum of DFT terms
CStilde_low = -inf(1,mean_out.mmax -l_star+1); %initialize various sums of DFT terms for necessary conditions
CStilde_up = inf(1,mean_out.mmax -l_star+1); %initialize various sums of DFT terms for necessary conditions
errest=zeros(mean_out.mmax -mean_out.mmin+1,1); %initialize error estimates
appxinteg=zeros(mean_out.mmax -mean_out.mmin+1,1); %initialize approximations to integral
exit_len = 2;
exit=false(1,exit_len); %we start the algorithm with all warning flags down

%% Initial points and FFT
mean_out.nSample=2^mean_out.mmin; %total number of points to start with
n0=mean_out.nSample; %initial number of points
xpts=mod(bsxfun(@plus, gail.lattice_gen(1,n0,mean_out.d), mean_out.shiftVal),1); %grab Lattice points

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
    
    if mean_out.nCV 
        evenval=yg(ptind, (1:mean_out.nCV ));
        oddval=yg(~ptind, (1:mean_out.nCV ));
        yg(ptind, (1:mean_out.nCV ))=(evenval+coefv.*oddval)/2;
        yg(~ptind, (1:mean_out.nCV ))=(evenval-coefv.*oddval)/2;
        
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
if mean_out.nCV 
    C=[ones(mean_out.nf,1); zeros(mean_out.nCV,1)];
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
mean_out.errBd=mean_out.inflate(mean_out.mmin)*Stilde(1);
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
if mean_out.nCV 
    q = mean(yval) - mu*beta(mean_out.nf+1:end,:);
else
    q =mean(yval);
end

% Check the end of the algorithm
q = q - errest(1)*(max(mean_out.absTol, mean_out.relTol*abs(q + errest(1)))...
    - max(mean_out.absTol, mean_out.relTol*abs(q - errest(1))))/...
    (max(mean_out.absTol, mean_out.relTol*abs(q + errest(1)))...
    + max(mean_out.absTol, mean_out.relTol*abs(q - errest(1)))); % Optimal estimator

q=q(1);
appxinteg(1)=q;

is_done = false;
if 4*errest(1)^2/(max(mean_out.absTol, mean_out.relTol*abs(q + errest(1)))...
        + max(mean_out.absTol, mean_out.relTol*abs(q - errest(1))))^2 <= 1
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
    xnext=mod(bsxfun(@plus, gail.lattice_gen(nnext+1,2*nnext,mean_out.d), mean_out.shiftVal),1);
    n0=n0+nnext;
    
    % check for using control variates or not
    if mean_out.nCV  == 0
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
            evenval=yg(ptind, (1:mean_out.nCV ));
            oddval=yg(~ptind, (1:mean_out.nCV ));
            yg(ptind, (1:mean_out.nCV ))=(evenval+coefv.*oddval)/2;
            yg(~ptind, (1:mean_out.nCV ))=(evenval-coefv.*oddval)/2;
        end
        
        X = yg(kappanumap(2^(m-r_lag-1)+1:end), (1:mean_out.nCV ));
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
    meff=m-mean_out.mmin+1;
    Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
    mean_out.errBd=mean_out.inflate(m)*Stilde(meff);
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
    if mean_out.nCV 
        q= mean(yval) - mu*beta(mean_out.nf+1:end,:);
    else
        q=mean(yval);
    end
    
    % Check the end of the algorithm
    q = q - errest(meff)*(max(mean_out.absTol, mean_out.relTol*abs(q + errest(meff)))...
        - max(mean_out.absTol, mean_out.relTol*abs(q - errest(meff))))/...
        (max(mean_out.absTol, mean_out.relTol*abs(q + errest(meff)))...
        + max(mean_out.absTol, mean_out.relTol*abs(q - errest(meff)))); % Optimal estimator
    
    appxinteg(meff)=q;
    
    if 4*errest(meff)^2/(max(mean_out.absTol, mean_out.relTol*abs(q + errest(meff)))...
            + max(mean_out.absTol, mean_out.relTol*abs(q - errest(meff))))^2 <= 1
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
mean_out.mu=q;

meanf=mean_out.mu;

end


