function [mu, out_param]=dt_meanMC_g_TrafficModel(tol)
% test traffic model
%
% dt_meanMC_g_TrafficModel(tol)
%    returns (mu)
%
% Examples:
%
% >> dt_meanMC_g_TrafficModel(1e-2)
% 
% ans =
% 
%     4.4***
%
% >> dt_meanMC_g_TrafficModel(1e-3)
%
% Warning: tried to evaluate at *** samples, which is more than
% the allowed maximum of *** samples. Just use the maximum sample
% budget.
% 
% >> dt_meanMC_g_TrafficModel('hi')
% ??? Error using ***
% 
in_param.tbudget=30;% given time budget 
in_param.nbudget=1e6;% given sample budget
in_param.n_sigma=70;% given intitial sample size to estimate sigma
in_param.fudge=1.1;%given fudge factor
in_param.alpha=0.01;% given uncertainty
in_param.abstol=tol;%given error tolerance
in_param.reltol = 0; % set relative tolerance 0
[mu, out_param] = meanMC_g(@Ytrafficmodel,in_param);
end
