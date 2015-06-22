%This script will construct the flat function.
% 
%  Generates Figure 4.1 in the thesis with a=0.5 and b=1
%
%  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for Univariate
%  Function Minimization. MS thesis, Illinois Institute of Technology,
%  2014.

function y=UniFunMin_Plot_Flat(a,b)

f=@(x) exp(-b./(x-a).^2);

%% Plot
[~,~,MATLABVERSION]=GAILstart(0);
if usejava('jvm') || MATLABVERSION <= 7.12
  t=0:0.000001:1;
  y=f(t);
  plot(t,f(t),'LineWidth',2)
  xlim([0 1])
  title('Example of a flat function')
  xlabel('X');ylabel('Y');
  
  gail.save_eps('UniFunMinOutput', 'UniFunMinPlotFlat');
end