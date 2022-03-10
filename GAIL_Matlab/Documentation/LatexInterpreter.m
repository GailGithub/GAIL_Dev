function LatexInterpreter
% LATEXINTERPRETER sets the interpreter for plot axes, text annotations,
% and labels to be LaTeX rather than TeX

%% From https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex?s_tid=mlc_ans_email_view#answer_803656

list_factory = fieldnames(get(groot,'factory')); 
index_interpreter = find(contains(list_factory,'Interpreter'));

for i = 1:length(index_interpreter) 
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default '); 
    set(groot, default_name,'latex '); 
end

%% Below is the old version
% if str2double(getfield(ver('MATLAB'),'Version')) >= 8.5 %the next part only works for later versions of MATLAB
%    set(groot,'defaultAxesTickLabelInterpreter','latex', ... %LaTeX interpreted axis tics
%       'defaultLegendInterpreter','latex', ... %LaTeX interpreted legends
%       'defaultTextArrowShapeInterpreter','latex', ... %LaTeX interpreted annotations
%       'defaultTextBoxShapeInterpreter','latex', ... %LaTeX interpreted annotations
%       'defaultColorbarTickLabelInterpreter','latex')  %LaTeX interpreted colorbar
%    if isprop(groot,'defaultPolaraxesTickLabelInterpreter')
%       set(groot,'defaultPolaraxesTickLabelInterpreter','latex')
%    end
% end
% set(groot,'defaultTextInterpreter','latex') %LaTeX interpreted axis labels

