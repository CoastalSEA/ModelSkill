ax = gca;
xmnmx = xlim(ax);
ymnmx = ylim(ax);
hold on

y = 0.0004*(xmnmx).^0.91;
plot(ax,xmnmx,y,'--k','DisplayName','O''Brien-Jarrett')

y = 0.0123*(xmnmx).^0.81;
plot(ax,xmnmx,y,'-.g','DisplayName','T-H-H Rias')

y = 0.0722*(xmnmx).^0.81;
plot(ax,xmnmx,y,':b','DisplayName','T-H-H Fjords')

ylim(ymnmx);

hold off