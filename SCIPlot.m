function h = SCIPlot()
    grid off;
    box on;
    set(gca,'Linewidth',2);
    a=gca;
    a.FontSize=16;
    a.FontName='Helvetica';
    a.Box='on';
    set(gca,'looseInset',[0.01 0.01 0.01 0.01])
    set(gcf,'color','w');
    set(gca,'linewidth',1.5);
end

