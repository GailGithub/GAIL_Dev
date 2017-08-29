function [mu, out_param]=dt_meanMC_g_TrafficModel(tol)
% test traffic model
%
% dt_meanMC_g_TrafficModel(tol)
%    returns (mu)
%
% Examples:
%
% >> dt_meanMC_g_TrafficModel('hi')
% ??? Error using ***
%
% >> dt_meanMC_g_TrafficModel(1e-2)
% 
% ans =
% 
%     4.4***
%
% >> dt_meanMC_g_TrafficModel(1e-3)
%
% Warning: In order to achieve the guaranteed accuracy, at step *** , tried to
% evaluate at *** samples, which is more than the remaining ***
% samples. We will use all the samples left to estimate the mean without
% guarantee. ***
%
% 
in_param.tbudget=30;% given the time budget 
in_param.nbudget=1e6;% given the sample budget
in_param.nSig=70;% given the intitial sample size to estimate sigma
in_param.n1 = 70;% the initial sample size to estimate the mean
in_param.fudge=1.1;%given the fudge factor
in_param.alpha=0.05;% given the uncertainty 
in_param.abstol=tol;%given error tolerance
in_param.reltol = 0; 
% set relative error tolerance zero, which means only require the error
% meets absolute error tolerance
[mu, out_param] = meanMC_g(@Ytrafficmodel,in_param);%Call Ytrafficmodel to get the results
end
