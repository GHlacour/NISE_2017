d=load('2D.par.dat');
w1=d(:,2);
w2=reshape(w1,sqrt(length(w1)),sqrt(length(w1)));
w=w2(:,1);
t=d(:,4);
tmp=reshape(t,length(w),length(w));

clevel=-9:.1:.9;
cmap=colormap2d(5);
figure(1);clf;
tmp=tmp/max(max(abs(tmp)));
contourf(w,w,(tmp),clevel,'LineWidth',.5);axis tight; axis square; colormap(cmap);
set (gca, 'TickLength', [0.035 0.035], 'TickDir', 'out');
set (gca, 'xtick',[1550 1600 1650 1700 1750],'LineWidth',.75,'FontSize',20);
set (gca, 'ytick',[1550 1600 1650 1700 1750],'LineWidth',.75,'FontSize',20);
axis([1550 1750 1550 1750])
hold on;
plot(w,w,'k');

