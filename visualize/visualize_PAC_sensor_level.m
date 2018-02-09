function visualize_PAC_sensor_level(file1)
%%
p2df = '/home/rouxf/prj/TC/matFilesRev/';
load([p2df,file1.name]);

%%
Lid = ft_channelselection({'ML*'},raw1.label);
Rid = ft_channelselection({'MR*'},raw1.label);
Zid = ft_channelselection({'MZ*'},raw1.label);

sel_idxL  = zeros(length(Lid),1);
for it = 1:length(Lid)
    sel_idxL(it) = find(strcmp(raw1.label,Lid(it)));
end;

sel_idxR  = zeros(length(Rid),1);
for it = 1:length(Rid)
    sel_idxR(it) = find(strcmp(raw1.label,Rid(it)));
end;

sel_idxZ  = zeros(length(Zid),1);
for it = 1:length(Zid)
    sel_idxZ(it) = find(strcmp(raw1.label,Zid(it)));
end;
%%
[i1,i2,i3] = ind2sub(size(MI(sel_idxL,:,:)),find(MI(sel_idxL,:,:) == max(max(max(MI(sel_idxL,:,:))))));
[i4,i5,i6] = ind2sub(size(MI(sel_idxR,:,:)),find(MI(sel_idxR,:,:) == max(max(max(MI(sel_idxR,:,:))))));
[i7,i8,i9] = ind2sub(size(MI(sel_idxZ,:,:)),find(MI(sel_idxZ,:,:) == max(max(max(MI(sel_idxZ,:,:))))));

i1 = sel_idxL(i1);
i4 = sel_idxR(i4);
i7 = sel_idxZ(i7);

%%
figure;

subplot(1,3,[1]);
imagesc(pfoi,afoi,squeeze(mean(MI,1))');
axis xy;
ca = caxis;
xlabel('Frequency [Hz]','Fontsize',14);
ylabel('Frequency [Hz]','Fontsize',14);
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','MI [a.u.]','Fontsize',14);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
set(cb,'YTick',[min(get(cb,'YTick')) max(get(cb,'YTick'))]);

set(gca,'YTick',[afoi(1):40:afoi(length(afoi))]);
set(gca,'XTick',[pfoi(1):8:pfoi(length(pfoi))]);
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'ticklength',get(gca,'ticklength')*2);
set(gca,'FontName','Arial');

subplot(132);
hold on;
plot(pbins,squeeze(PAH(i1,i2,i3,:)),'b-','LineWidth',3);
plot(pbins,squeeze(PAH(i4,i5,i6,:)),'r-','LineWidth',3);
axis tight;
xlabel('\alpha-phase [rad]','Fontsize',14);
ylabel('\gamma-power [a.u.]','Fontsize',14);
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'FontName','Arial');
set(gca,'ticklength',get(gca,'ticklength')*2);
set(gca,'XTick',[-pi -pi/2 0 pi/2 pi]);
set(gca,'XTickLabel',{'-\pi' '-\pi/2' '0' '\pi/2' '\pi'},'Fontsize',14);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);

subplot(1,3,3);
dum = [];
dum.avg = squeeze(MI(:,i2,i3));
dum.time = 1;
dum.label = raw1.label;
dum.dimord = 'chan_time';

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.parameter = 'avg';
cfg.comment = 'no';
cfg.marker = 'off';

ft_topoplotER(cfg,dum);
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','MI [a.u.]');
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
set(cb,'YTick',[min(get(cb,'YTick')) max(get(cb,'YTick'))]);

set(gcf,'Color','w');

%%
savepath = '/home/rouxf/prj/TC/figures/';
for it = 1
    pos = get(gcf,'PaperPosition');
    dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
    set(gcf,'PaperPositionMode','auto');
    if dim(2)>dim(1)
        set(gcf,'PaperOrientation','landscape');
    end;
    it
    print(it,'-dsvg','-opengl','-r600',[savepath,'MIatsensorLevel',num2str(it),'.svg']);
    print(it,'-dpdf','-opengl','-r600',[savepath,'MIatsensorLevel',num2str(it),'.pdf']);
    %plot2svg([savepath,'sourceModel_alpha_generator',num2str(it),'.svg'],it,'png');
end;
