% RunMVNPdeltaplot produces Figure 4.5 in Lan Jiang's thesis
load TestcubMCon-MVNP-N500d3abstol1e-05rel1e-05-2015-09-12-00-38-36.mat
%set(0,'defaultlinelinewidth',4)
deltavec = [0.01,0.1,1,10];
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
figure
boxplot(res.iidtimetot,deltavec,'datalim',[0 200],'symbol','k+','colors','k')
xlabel('Delta')
ylabel('Time(in seconds)')
set(gca,'Yscale','log','Ylim',[0.01 200])
gail.save_eps('LanThesisOutput','boxplotofMVNP');
% figure
% boxplot(log10(res.iidtimetot),deltavec)
% xlabel('Delta')
% ylabel('Time')
% ticval = [0.01 0.1 1 10 100]
% %ticval=[0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12 10.24 20.48 60.96 121.92];
% set(gca,'YTick',log10(ticval),'YTickLabel',ticval,'Ylim',[-2 2.2])