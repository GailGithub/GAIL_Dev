function plotTestcubMCblack(plotTest,param)
[~,~,~,MATLABVERSION] = GAILstart(false); 
if usejava('jvm') || MATLABVERSION <= 7.12
figure
ntot=length(plotTest.err);
standard=true(size(plotTest.err));
if isfield(plotTest,'kurtvec') && isfield(plotTest,'kurtmax')
    smallkurt=plotTest.kurtvec<=plotTest.kurtmax;
    standard=standard &~smallkurt;
    scatter(plotTest.err(smallkurt,:),plotTest.time(smallkurt,:),plotTest.ptsize,'k*', ...
        'linewidth',1);
    hold on;
end

if isfield(plotTest,'exit')
   bigsample=plotTest.exit==1;
   standard=standard&~bigsample;
   scatter(plotTest.err(bigsample,:),plotTest.time(bigsample,:),plotTest.ptsize,'kd','filled');
   
end
erraug=[plotTest.errlowlimit; sort(plotTest.err); plotTest.errhilimit];
timeaug=[plotTest.timelowlimit; sort(plotTest.time); plotTest.timehilimit];
probaug=[0; ((0:ntot-1)'+1/2)/ntot; 1];
scatter(plotTest.err(standard,:),plotTest.time(standard,:),plotTest.ptsize,'k.')
ax1 = gca;
  
set(ax1,'XColor','k','YColor','k',...
    'XLim',[plotTest.errlowlimit plotTest.errhilimit],...
    'YLim',[plotTest.timelowlimit plotTest.timehilimit],...
    'LineWidth',plotTest.linewidth, ...
    'XScale','log','Yscale','log', ...
    'XTick',10.^(plotTest.logerrlo:plotTest.logerrhi))
    %'XGrid','on','YGrid','on','GridLineStyle','--',...

line([param.tol param.tol],[plotTest.timelowlimit plotTest.timehilimit],...
    'color','k','linestyle','--','linewidth',plotTest.linewidth)
% wherexsuccess=-0.02+0.75*(log10(param.tol)-plotTest.logerrlo)/...
%     (plotTest.logerrhi-plotTest.logerrlo);
% wherexfailure=0.28+0.75*(log10(param.tol)-plotTest.logerrlo)/...
%     (plotTest.logerrhi-plotTest.logerrlo);
% whereysuccess=0.2;
% whereyfailure=0.85;
% annotation('textarrow',[wherexsuccess,0.15],whereysuccess*[1 1],...
%     'String','Success','FontWeight','bold',...
%     'linewidth',2);
% annotation('textarrow',[wherexfailure,0.88],whereyfailure*[1 1],...
%     'String','Failure','FontWeight','bold',...
%     'linewidth',2);
xlabel('Error')
ylabel('Time (seconds)')
axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','w','YColor','k',...
           'XLim',[plotTest.errlowlimit plotTest.errhilimit],...
           'Xscale','log','Xtick',[],...
           'YLim',[0 1],'Linewidth',plotTest.linewidth);
%ylabel('Probability')
line(erraug,probaug,'color','k','linestyle','-','linewidth',plotTest.linewidth)
axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','k','YColor','k',...
           'YLim',[plotTest.timelowlimit plotTest.timehilimit],...
           'Yscale','log','Ytick',[],...
           'XLim',[0 1],'Linewidth',plotTest.linewidth);
%xlabel('Probability')
line(probaug,timeaug,'color','k','linestyle','-.','linewidth',plotTest.linewidth)
gail.save_eps('MCQMC2012PaperOutput/Results',plotTest.name);
end
end