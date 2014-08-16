function eps = tolfun(abstol,reltol,theta,mu,errtype)
% TOLFUN generalized error tolerance function.
%
% Input Parameters:
% abstol --- absolute error tolertance
% reltol --- relative error tolerance
% theta --- parameter in 'theta' case
% mu --- true mean
% errtype --- different option of tolerance function

switch errtype
    case 'theta' % the theta case
        %theta=0---absolute error
        %theta=1---relative error
        eps  = theta*abstol+(1-theta)*reltol;
    case 'max' % the max case
        eps  = max(abstol,reltol*abs(mu));
end
end