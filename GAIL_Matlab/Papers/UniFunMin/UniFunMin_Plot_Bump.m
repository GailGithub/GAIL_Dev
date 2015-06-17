%This function will construct the example of the bump test functions.
% 
%  Generates Figure 3.1 in the thesis with a=0.03 and z=0.4
%
%  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for Univariate
%  Function Minimization. MS thesis, Illinois Institute of Technology,
%  2014.


function y=UniFunMin_Plot_Bump(a,z)

f=@(x) 0.5/a^2*(-4*a^2-(x-z).^2-(x-z-a).*abs(x-z-a)+(x-z+a).*...
    abs(x-z+a)).*(x>=z-2*a).*(x<=z+2*a);

%% Plot
[~,~,MATLABVERSION]=GAILstart(0);
[~,~,MATLABVERSION]=GAILstart(0);
if usejava('jvm') || MATLABVERSION <= 7.12
  t=0:0.00001:1;
  y=f(t);
  plot(t,y,'LineWidth',2)
  xlim([0 1])
  ylim([-1.2 0.2])
  title('An example of the bump test functions')
  xlabel('X');ylabel('Y');
  
  
  %% Save output
  gail.save_eps('UniFunMinOutput', 'UniFunMinPlotBump');
end