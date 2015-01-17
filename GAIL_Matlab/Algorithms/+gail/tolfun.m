function eps = tolfun(abstol,reltol,theta,mu,toltype)
% TOLFUN generalized error tolerance function.
%
% Input Parameters:
% abstol --- absolute error tolertance
% reltol --- relative error tolerance
% theta --- parameter in 'theta' case
% mu --- true mean
% toltype --- different options of tolerance function

switch toltype
    case 'comb' % the linear combination of two tolerances
        %theta=0---relative error tolarance
        %theta=1---absolute error tolerance
        eps  = theta*abstol+(1-theta)*reltol*abs(mu);
    case 'max' % the max case
        eps  = max(abstol,reltol*abs(mu));
end
end