function [muhat,out]=cubMLELattice(f,d,absTol,relTol,order,ptransform)
%CUBMLE Monte Carlo method to estimate the mean of a random variable
%
%   tmu = cubMLELattice(f,absTol,relTol,alpha,nSig,inflate) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples may be of one of several kinds.  The default values are n=2^10 and
%   d = 1 Input f is a function handle that accepts an n x d matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%
%
% This is a heuristic algorithm based on a Central Limit Theorem
% approximation
if nargin < 6
    ptransform = 'Baker';
    if nargin < 5
        order = 2; %type of sampling, scrambled Sobol
        if nargin < 4
            relTol = 0;
            if nargin < 3
                absTol = 0.01;
                if nargin < 2
                    d = 1; %dimension
                    if nargin < 1
                        f = @(x) x.^2; %function
                    end
                end
            end
        end
    end
end

tstart = tic; %start the clock
z = [1, 364981, 245389, 97823, 488939, 62609, 400749, 385317, 21281, 223487]; % generator from Hickernell's paper
%z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599]; %generator
z = z(1:d);
mmin = 10;
mmax = 23;
mvec = mmin:mmax;
numM = length(mvec);
errorBdAll = zeros(length(mvec),1);
muhatAll = zeros(length(mvec),1);
aMLEAll = zeros(length(mvec),1);
shift = rand(1,d);
% ff = f; %no transformation
% ff = @(x) f(1 - abs(1-2*x)); %folding transformation
ff = PeriodTx(f, ptransform);

for ii = 1:numM
    m = mvec(ii);
    n = 2^m
    
    %Update function values
    if ii == 1
        brIndices = bitrevorder((0:1/n:1-1/n));
        xun = mod(bsxfun(@times,brIndices',z),1);  % unshifted
        x = mod(bsxfun(@plus,xun,shift),1);  % shifted
        fx = ff(x);
        ftilde = fft(bitrevorder(fx))/n;
        
        %% Efficient FFT computation algorithm
        ftildeNew=ff(x); %evaluate integrand
        
        %% Compute initial FFT
        for l=0:mmin-1
            nl=2^l;
            nmminlm1=2^(mmin-l-1);
            ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
            coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
            coefv=repmat(coef,nmminlm1,1);
            evenval=ftildeNew(ptind);
            oddval=ftildeNew(~ptind);
            ftildeNew(ptind)=(evenval+coefv.*oddval)/2;
            ftildeNew(~ptind)=(evenval-coefv.*oddval)/2;
        end
        
    else
        brIndices = bitrevorder((1/n:2/n:1-1/n));
        xunnew = mod(bsxfun(@times,brIndices',z),1);
        xnew = mod(bsxfun(@plus,xunnew,shift),1);
        if 0
            temp = zeros(n,d);
            temp(1:2:n-1,:) = xun;
            temp(2:2:n,:) = xunnew;
            xun = temp;    % save unshifted for further iterations
            temp(1:2:n-1,:) = x;
            temp(2:2:n,:) = xnew;
            x = temp; % saving the shifted
        else
            xun = [xun;xunnew];
            x = [x;xnew];
        end
        fnew = ff(xnew);
        %fx = reshape([fx fnew]',n,1);
        fx = [fx;fnew];
        ftilde = fft(bitrevorder(fx))/n;
        
        %% Efficient FFT computation algorithm
        mnext=m-1;
        ftildeNextNew=ff(xnew);  % initialize for inplace computation
        
        %% Compute initial FFT on next points
        for l=0:mnext-1
            nl=2^l;
            nmminlm1=2^(mnext-l-1);
            ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
            coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
            coefv=repmat(coef,nmminlm1,1);
            evenval=ftildeNextNew(ptind);
            oddval=ftildeNextNew(~ptind);
            ftildeNextNew(ptind)=(evenval+coefv.*oddval)/2;
            ftildeNextNew(~ptind)=(evenval-coefv.*oddval)/2;
        end
        
        %% Compute FFT on all points
        ftildeNew=[ftildeNew;ftildeNextNew];
        nl=2^mnext;
        ptind=[true(nl,1); false(nl,1)];
        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
        coefv=repmat(coef,nmminlm1,1);
        evenval=ftildeNew(ptind);
        oddval=ftildeNew(~ptind);
        ftildeNew(ptind)=(evenval+coefv.*oddval)/2;
        ftildeNew(~ptind)=(evenval-coefv.*oddval)/2;
    end
    if sum((abs(ftildeNew-ftilde)))/n > 0.000001
        fprintf('FFT values differ too much')
    end
    ftilde = ftildeNew;
    
    %% figure(1); loglog(abs(ftilde)); figure(2); loglog((abs(ftildeNew)))
    
    br_xun = (xun);
    %Compute MLE parameter
    lnaMLE = fminbnd(@(lna) ...
        MLEKernel(exp(lna),br_xun,ftilde,order), ...
        -5,0,optimset('TolX',1e-2)); % -5,5
    aMLE = exp(lnaMLE);
    [loss,Ktilde,RKHSnormSq,K] = MLEKernel(aMLE,br_xun,ftilde,order);
    wt = 1./Ktilde(1);
    
    %Check error criterion
    DSC = abs(1 - wt); %/n;
    %DSC = 1 - n/sum(K);
    %DSC = (1/n) - 1/sum(Ktilde); % wrong
    out.ErrBd = 2.58*sqrt(DSC)*RKHSnormSq;
    muhat = ftilde(1)/Ktilde(1);
    muminus = muhat - out.ErrBd;
    muplus = muhat + out.ErrBd;
    muhatAll(ii) = muhat;
    errorBdAll(ii) = out.ErrBd;
    aMLEAll(ii) = aMLE;
    
    if 2*out.ErrBd <= ...
            max(absTol,relTol*abs(muminus)) + max(absTol,relTol*abs(muplus))
        
        fprintf('%d Error bound met %e !\n', n, out.ErrBd)
        if errorBdAll(ii)==0
            errorBdAll(ii) = eps;
        end
        % break
    end
    
end
out.n = n;
out.time = toc(tstart);
out.ErrBdAll = errorBdAll;
out.muhatAll = muhatAll;
out.mvec = mvec;
out.aMLEAll = aMLEAll;
end

function [loss,Ktilde,RKHSnormSq,K] = MLEKernel(a,xun,ftilde,order)

constMult = -(-1)^(order/2)*(2*pi)^order/factorial(order);
if order == 2
    K = prod(1 + a*constMult*(-xun.*(1-xun) + 1/6),2);
elseif order == 4
    K = prod(1 + a*constMult*(xun.^2.*(xun.*(xun-2) +1) - 1/30),2);
else
    error('Bernoulli order not implemented !');
end
n =  length(K);
Ktilde = real(fft_DIT((K)));
%Ktilde = real(fft_DIT(K));


%Ktilde = real(fft(K-1))/n; Ktilde(1) = Ktilde(1) + 1; % this is done to improve accuracy, to reduce zero values


%RKHSnorm = mean(abs(ftilde).^2./Ktilde);
%RKHSnorm = mean(abs(ftilde(Ktilde~=0)).^2./Ktilde(Ktilde~=0));
RKHSnormSq = mean(abs(ftilde(Ktilde~=0))./sqrt(Ktilde(Ktilde~=0)));


loss = mean(log(Ktilde(Ktilde~=0))) + 2*log(RKHSnormSq);
a;
ftilde(1)/Ktilde(1);
end

function f = PeriodTx(fInput, ptransform)

if strcmp(ptransform,'Baker')
    f=@(x) fInput(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(ptransform,'C0')
    f=@(x) fInput(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(ptransform,'C1')
    f=@(x) fInput(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(ptransform,'C1sin')
    f=@(x) fInput(x-sin(2*pi*x)/(2*pi)).*prod(2*sin(pi*x).^2,2); % Sidi C^1 transform
elseif strcmp(ptransform,'C2sin')
    psi3 = @(t) (8-9*cos(pi*t)+cos(3*pi*t))/16; psi3_1 = @(t) (9*sin(pi*t)*pi- sin(3*pi*t)*3*pi)/16;
    f=@(x) fInput(psi3(x)).*prod(psi3_1(x),2);
elseif strcmp(ptransform,'C3sin')
    psi4 = @(t) (12*pi*t-8*sin(2*pi*t)+sin(4*pi*t))/(12*pi); psi4_1 = @(t) (12*pi-8*cos(2*pi*t)*2*pi+sin(4*pi*t)*4*pi)/(12*pi);
    f=@(x) fInput(psi4(x)).*prod(psi4_1(x),2);
elseif strcmp(ptransform,'none')
    % do nothing
    f=@(x) fInput(x);
else
    error('Error: Periodization transform %s not implemented', ptransform);
end
end


function y = fft_DIT( y )
nmmin = log2(length(y));
%y = bitrevorder(y);
for l=0:nmmin-1
    nl=2^l;
    nmminlm1=2^(nmmin-l-1);
    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
    coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
    coefv=repmat(coef,nmminlm1,1);
    evenval=y(ptind);
    oddval=y(~ptind);
    y(ptind)=(evenval+coefv.*oddval)/2;
    y(~ptind)=(evenval-coefv.*oddval)/2;
end
end
