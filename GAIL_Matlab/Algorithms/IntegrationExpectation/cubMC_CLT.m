function [q, out_param] = cubMC_CLT(varargin)
t_start = tic;

mean_inp = gail.cubMCParam(cell2mat(varargin)); %parse the input and check it for errors
mean_out = gail.cubYOut(mean_inp); %create the output class

[q, out_param]=meanMC_CLT(mean_out.Y, mean_out.err.absTol, ...
   mean_out.err.relTol, mean_out.alpha, mean_out.nSig, mean_out.CM.inflate)

end 