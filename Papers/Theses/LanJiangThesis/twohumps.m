% twohumps is the function used to plot Figure 1.1 in Lan Jiang's thesis
% PlotTwoHumps.m calls this function

function f = twohumps(x)
f = 1/(sqrt(2*pi))*(exp(-x.^2/2)*0.99+exp(-(x-200).^2/2)*0.01);
end


