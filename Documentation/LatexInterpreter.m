function LatexInterpreter
% LATEXINTERPRETER sets the interpreter for plot axes, text annotations,
% and labels to be LaTeX rather than TeX
if str2double(getfield(ver('MATLAB'),'Version')) >= 8.5 %the next part only works for later versions of MATLAB
   set(groot,'defaultAxesTickLabelInterpreter','latex', ... %LaTeX interpreted axis tics
      'defaultLegendInterpreter','latex', ... %LaTeX interpreted legends
      'defaultTextArrowShapeInterpreter','latex', ... %LaTeX interpreted annotations
      'defaultTextBoxShapeInterpreter','latex', ... %LaTeX interpreted annotations
      'defaultColorbarTickLabelInterpreter','latex')  %LaTeX interpreted colorbar
   if isprop(groot,'defaultPolaraxesTickLabelInterpreter')
      set(groot,'defaultPolaraxesTickLabelInterpreter','latex')
   end
end
set(groot,'defaultTextInterpreter','latex') %LaTeX interpreted axis labels

