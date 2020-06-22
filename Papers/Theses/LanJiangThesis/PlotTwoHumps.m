% This script plots Figure 1.1 in Lan Jiang's thesis
% Function called: twohumps

A = -20:0.01:250;
B=twohumps(A);
plot(A,B,'linewidth',2)
%ylim([10^(-15) 1])  
xlabel('y')
ylabel('Probability density function of y')
set(gca,'Yscale','log','Ylim',[10^(-15) 1])
gail.save_eps('LanThesisOutput/mixturegaussian','mixturegaussianpdf');