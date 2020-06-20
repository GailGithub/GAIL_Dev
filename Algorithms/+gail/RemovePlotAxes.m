%% Remove all axes in a plot
hFig = figure;
color = get(hFig,'Color');
set(gca,'XColor',color,'YColor',color,'TickDir','out')
h = plot([1 2],'k');
delete(h)
set(gca,'Visible','off')
hold on

