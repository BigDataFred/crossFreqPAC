%%
figure;
subplot(321);
a = gca;
hold on;
%plot(dleM(:,1),'bs-','LineWidth',3,'MarkerFaceColor','b');
%plot(dleM(:,2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleM(:,1:2),2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleM2(:,1:2),2),'ks-','LineWidth',3,'MarkerFaceColor','k');
plot(mean(dleM3(:,1:2),2),'bs-','LineWidth',3,'MarkerFaceColor','b');

title('Parietal cortex');
ylabel('DLEm [mm]');
axis tight;
set(gca,'xTick',[1:size(dleM,1)]);
set(gca,'Xdir','reverse');
set(gca,'XTickLabel',SNR(1:4));
xlabel('SNR ');

subplot(322);
a = [a gca];
hold on;
%plot((dleM(:,3)),'bs-','LineWidth',3,'MarkerFaceColor','b');
%plot((dleM(:,4)),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleM(:,3:4),2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleM2(:,3:4),2),'ks-','LineWidth',3,'MarkerFaceColor','k');
plot(mean(dleM3(:,3:4),2),'bs-','LineWidth',3,'MarkerFaceColor','b');
title('Thalamus');
ylabel('DLEm [mm]');
axis tight;
set(gca,'xTick',[1:size(dleM,1)]);
set(gca,'Xdir','reverse');
set(gca,'XTickLabel',SNR(1:size(dleM,1)));
xlabel('SNR');

subplot(323);
a = [a gca];
hold on;
%plot(dleG(:,1),'bs-','LineWidth',3,'MarkerFaceColor','b');
%plot(dleG(:,2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleG(:,1:2),2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleG2(:,1:2),2),'ks-','LineWidth',3,'MarkerFaceColor','k');
plot(mean(dleG3(:,1:2),2),'bs-','LineWidth',3,'MarkerFaceColor','b');
ylabel('DLEg [mm]');
axis tight;
set(gca,'xTick',[1:size(dleM,1)]);
set(gca,'Xdir','reverse');
set(gca,'XTickLabel',SNR(1:size(dleM,1)));
xlabel('SNR');

subplot(324);
a = [a gca];
hold on;
%plot((dleG(:,3)),'bs-','LineWidth',3,'MarkerFaceColor','b');
%plot((dleG(:,4)),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleG(:,3:4),2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(dleG2(:,3:4),2),'ks-','LineWidth',3,'MarkerFaceColor','k');
plot(mean(dleG3(:,3:4),2),'bs-','LineWidth',3,'MarkerFaceColor','b');
ylabel('DLEg [mm]');
axis tight;
set(gca,'xTick',[1:size(dleM,1)]);
set(gca,'Xdir','reverse');
set(gca,'XTickLabel',SNR(1:size(dleM,1)));
xlabel('SNR');

subplot(325);
a = [a gca];
hold on;
plot(mean(AUC(:,1:2),2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(AUC2(:,1:2),2),'ks-','LineWidth',3,'MarkerFaceColor','k');
plot(mean(AUC3(:,1:2),2),'bs-','LineWidth',3,'MarkerFaceColor','b');
ylabel('AUC [a.u.]');
axis tight;
set(gca,'xTick',[1:size(AUC,1)]);
set(gca,'Xdir','reverse');
set(gca,'XTickLabel',SNR(1:size(AUC,1)));
xlabel('SNR');

subplot(326);
a = [a gca];
hold on;
plot(mean(AUC(:,3:4),2),'rs-','LineWidth',3,'MarkerFaceColor','r');
plot(mean(AUC2(:,3:4),2),'ks-','LineWidth',3,'MarkerFaceColor','k');
plot(mean(AUC3(:,3:4),2),'bs-','LineWidth',3,'MarkerFaceColor','b');
ylabel('AUC [a.u.]');
axis tight;
set(gca,'xTick',[1:size(AUC,1)]);
set(gca,'Xdir','reverse');
set(gca,'XTickLabel',SNR(1:size(AUC,1)));
xlabel('SNR');

for it = 1:length(a)
    set(get(a(it),'YLabel'),'Fontsize',14);
    if it <5
        set(a(it),'YLim',[0 max(get(a(it),'YLim'))]);
    else
        set(a(it),'YLim',[0 1]);
    end;
    try set(a(it),'YTick',[5 floor(max(get(a(it),'YLim')))]);catch;set(a(it),'YTick',[0 1]);end;
end;
set(a,'LineWidth',3);
set(a,'Fontsize',12);
set(a,'FontName','Arial');
set(a,'ticklength',unique(cell2mat(get(a,'ticklength')))'*2);
set(gcf,'Color','w');

%%
savepath = '/home/rouxf/prj/TC/figures/';

pos = get(gcf,'PaperPosition');
dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
set(gcf,'PaperPositionMode','auto');
if dim(2)>dim(1)
    set(gcf,'PaperOrientation','landscape');
end;

print(gcf,'-dsvg','-opengl','-r1000',[savepath,'DLEsourceLevel',num2str(it),'.svg']);
print(gcf,'-dpdf','-opengl','-r1000',[savepath,'DLEsourceLevel',num2str(it),'.pdf']);
