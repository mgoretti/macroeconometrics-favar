function formatGraph(xlab, ylab)
% The function saveGraph(path, hxlabel, hylabel, ratio) makes the plot look 
% nice and increase font size etc.

hx = xlabel(xlab);
hy = ylabel(ylab);
xlim([1 20])
grid on;
set(gca,'fontsize',14,'fontname','Helvetica','box','off','tickdir','out','ticklength',[.02 .02],'xcolor',0.5*[1 1 1],'ycolor',0.5*[1 1 1]);
set([hx; hy],'fontsize',12,'fontname','avantgarde','color',[.3 .3 .3]);

end

