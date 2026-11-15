function export_pdf(fig,path)
%EXPORT_PDF Vector PDF, separate from all scientific calculations.
set(fig,'Color','w');
axes=findall(fig,'Type','axes');
for k=1:numel(axes)
 set(axes(k),'FontName','Times New Roman','FontSize',11,'LineWidth',.8,'TickDir','out');
 set(axes(k),'XMinorGrid','off','YMinorGrid','off','GridAlpha',.15);
 box(axes(k),'off');
end
exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white');
end
