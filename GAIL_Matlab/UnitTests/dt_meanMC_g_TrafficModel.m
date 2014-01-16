function mu=dt_meanMC_g_TrafficModel(tol)
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
% >> dt_meanMC_g_TrafficModel(5e-3)
% 
% ans =
% 
%     4.45***
% 
% >> dt_meanMC_g_TrafficModel('hi')
% ??? Error using ***meanMC_g***Argument 'abstol' failed validation isnumeric.
% 
%
in_param.timebudget=15;% given time budget 
in_param.nbudget=1e6;% given sample budget
in_param.n0=70;% given intitial sample size to estimate sigma
in_param.fudge=1.1;%given fudge factor
in_param.alpha=0.01;% given uncertainty
in_param.abstol=tol;%given error tolerance
in_param.npcmax=1e6;% given piecewise maximum to calculate the mu
mu=meanMC_g(@Ytrafficmodel,in_param);
end
