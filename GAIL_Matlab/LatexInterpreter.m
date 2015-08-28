function LatexInterpreter
% LATEXINTERPRETER sets the interpreter for plot axes, text annotations,
% and labels to be LaTeX rather than TeX
if str2double(getfield(ver('MATLAB'), 'Version')) >= 8.5 %the next part only works for later versions of MATLAB
   set(0, 'defaultAxesTickLabelInterpreter','latex', ... %LaTeX interpreted axis tics
      'defaultLegendInterpreter','latex'); %LaTeX interpreted legends
end
set(0,'defaultTextInterpreter','latex') %LaTeX interpreted axis labels

