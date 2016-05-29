function saveGraph(path, hxlabel, hylabel, ratio)
% The function saveGraph(path, hxlabel, hylabel, ratio) makes the plot look 
% nice and increase font size etc.

hx = xlabel(hxlabel);
hy = ylabel(hylabel);

set(gca,'fontsize',14,'fontname','Helvetica','box','off','tickdir','out','ticklength',[.02 .02],'xcolor',0.5*[1 1 1],'ycolor',0.5*[1 1 1]);
set([hx; hy],'fontsize',12,'fontname','avantgarde','color',[.3 .3 .3]);
%grid on;

hold off;
w = 10*ratio; h = 10;
set(gcf, 'PaperPosition', [0 0 w h]); %Position plot at left hand corner with width w and height h.
set(gcf, 'PaperSize', [w h]); %Set the paper to have width w and height h.
saveas(gcf, path, 'pdf') %Save figure

end
