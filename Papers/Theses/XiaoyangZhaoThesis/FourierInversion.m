function y = FourierInversion( h,Vs,Vt,dt,nu,kappa,theta,N,x,start )
%FourierIncersion
%Matlab implementation of the FourierInversion function used to invert the characteristic
%function needed to apply the exact sampling numerically

% Code source:Kienitz, J., & Wetterau, D. (n.d.). Financial Modelling 
% Theory Implementation and Practice with MATLAB (pp. 296).
A = h*x/pi;
vec = (start:N)*h;
B = sum((2/pi)*(sin(vec.*x)./(start:N)).*...
    exp(CharFun(vec,Vs,Vt,dt,nu,kappa,theta)));
y = A+B;

end

