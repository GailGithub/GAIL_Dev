function [ Price ] = ...
    ExactSampling_Heston(S0,V0,K,r,T,kappa,theta,nu,rho,Nt,NSim)
%Matlab implementation of the Exact Sampling scheme of Broadie and Kaya for the Heston
%stochastic volatility model

% Code source:Kienitz, J., & Wetterau, D. (n.d.). Financial Modelling 
% Theory Implementation and Practice with MATLAB (pp. 296).

dt = T/Nt;                                         % equidistant step
t = 0:dt:T;                                        % time grid
N = 800;
h = 0.5;
start = 1;
V = ones(NSim,Nt+1);
V = V0*V(:,1);                                     % init variance
S = ones(NSim,Nt+1);                               
S = S0*S(:,1);                                     % init asset
Inverse = zeros(NSim,1);                           % for num inversion
options = optimset('TolX',1e-8);                   % for num inversion

for i = 2:Nt+1
    dt = t(i)-t(i-1);                              % time step
    A = 4*kappa*theta/(nu^2);                      % parameter chi^2
    B = (4*kappa*exp(-kappa*dt)*V(:,i-1))...
        /((nu^2)*(1-exp(-kappa*dt)));              % parameter chi^2
    C = (((nu^2)*(1-exp(-kappa*dt)))/(4*kappa));   % parameter chi^2
    V(:,i) = C*ncx2rnd(A,B,NSim,1);                % variate for V
    
    uniforms = rand(NSim,1);                       % uniforms
    
    % apply numerical inversion
    for j = 1:NSim
        Inverse(j) = fzero(@(x) FourierInversion(h,V(j),V(j,i),...
            dt,nu,kappa,theta,N,x,start)-uniforms(j),0.05,options);
    end
    % integral
    Integral = (1/nu)*(V(:,i)-V(:,i-1)-kappa*theta*(dt)+kappa*Inverse);
    mu = log(S(:,i-1))+r*dt-0.5*Inverse+rho*Integral;
    sig = (1-rho^2)*Inverse;
    S(:,i) = exp(mu+randn(NSim,1).*sqrt(sig));      % variate from path
end
Price = mean(max(S(:,end)-K,0))*exp(-r*T);
    
    
end

