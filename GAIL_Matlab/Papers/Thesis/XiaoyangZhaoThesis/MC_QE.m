function [pathS, pathV] = MC_QE(S0,r,d,T,Vinst,Vlong,kappa,nu,rho,...
    NTime,NSim,NBatches)
%************************************************************************
%
% Code source:Kienitz, J., & Wetterau, D. (n.d.). Financial Modelling 
% Theory Implementation and Practice with MATLAB (pp. 298-299). 
%
%************************************************************************
% discretization for the Heston model
% using QE scheme

dT = T/NTime;                           % time step

pathS = zeros(NSim,NTime+1,NBatches);   % output pathS
pathV = zeros(NSim,NTime+1,NBatches);   % output pathV

lnS1 = zeros(NSim,NTime+1);     % logspot price path
lnS1(:,1)=log(S0*exp(-d*T));    % set S(0) adjust with dividend

V2 = zeros(NSim,NTime+1);       % Variance path
V2(:,1) = Vinst;                % set V0

k1 = exp(-kappa*dT);
k2 = nu^2*k1.*(1-k1)/kappa;
k3 = exp(kappa*dT)*0.5.*k2.*(1-k1).*Vlong;

psiC = 1.5;                      % psiC in (1,2)
gamma1 = .5;                     % For PredictorCorrector
gamma2 = .5;                     % For PredictorCorrector

c1 = (r-d)*dT;                   %adjustment due to drift
c2 = -rho*kappa*Vlong*dT/nu;%used to determine K0

K0 = c1 + c2;                    % drift adjusted K0
K1 = gamma1*dT*(kappa*rho/nu - .5)-rho/nu; %K1
K2 = gamma2*dT*(kappa*rho/nu - .5)+rho/nu; %K2
K3 = gamma1*dT*(1-rho^2);                            %K3
K4 = gamma2*dT*(1-rho^2);                            %K4

for l = 1:NBatches
    UV1 = rand(NSim,NTime);     % uniforms
    dW2 = randn(NSim,NTime);    % Gaussians
    
    for i=2:NTime+1             % time loop
        m = Vlong + (V2(:,i-1)-Vlong)*k1;% mean (moment matching)
        s2 = V2(:,i-1)*k2 + k3;          % var (moment matching)
        psi = s2./m.^2;                  % psi compared to psiC
        
        psihat = 1./psi;
        b2 = 2*psihat - 1 + sqrt(2*psihat.*(2*psihat-1));
        a = m ./ (1 + b2);
        
        I1 = find(psi<=psiC); % Non-Central Chi^2 approx for psi < psiC
        if isempty(I1)
        else
            V2(I1,i) = a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
        end
        p = (psi - 1)./(psi + 1);               % for switching rule
        V2((UV1(:,i-1)<=p) & (psi>psiC),i) = 0; % case u<=p & psi>psiC
        I1b = find((UV1(:,i-1)>p) & (psi>psiC));% find is faster here!
        
        beta = (1 - p)./m;                      % for switching rule
        if isempty(I1b)
        else    % Psi^(-1)
            V2(I1b,i) = log((1-p(I1b))./(1-UV1(I1b,i-1)))./beta(I1b);
        end
        % log Euler Predictor-Corrector step
        lnS1(:,i) = lnS1(:,i-1) + K0 + K1.*V2(:,i-1) + K2.*V2(:,i)...
                  + sqrt(K3.*V2(:,i-1) + K4.*V2(:,i)).*dW2(:,i-1);
    end
    pathS(:,:,1) = exp(lnS1);
    pathV(:,:,1) = V2;
end
