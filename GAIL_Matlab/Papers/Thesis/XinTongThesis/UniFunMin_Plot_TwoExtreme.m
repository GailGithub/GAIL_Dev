%This function will construct the function with two local minimum points.
%  
%  Generates Figure 3.2 in the thesis with a1=0.3 and a2=0.75
%
%  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for Univariate
%  Function Minimization. MS thesis, Illinois Institute of Technology,
%  2014.

function y=UniFunMin_Plot_TwoExtreme(a1,a2)

f=@(x) -5*exp(-(10*(x-a1)).^2)-exp(-(10*(x-a2)).^2);

%% Plot
[~,~,MATLABVERSION]=GAILstart(0);
if usejava('jvm') || MATLABVERSION <= 7.12
  t=0:0.00001:1;
  y=f(t);
  plot(t,y,'LineWidth',2)
  xlim([0 1])
  ylim([-6 1])
  title('A function with two local minimum points')
  xlabel('X');ylabel('Y');
  
  gail.save_eps('UniFunMinOutput', 'UniFunMinTwoExtreme');
end