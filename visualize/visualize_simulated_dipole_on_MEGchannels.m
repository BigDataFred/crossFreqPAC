function visualize_simulated_dipole_on_MEGchannels(file1,file2)

%%
p2fd = '~rouxf/prj/TC/matFilesRev/';

file1 = dir([p2fd,file1]);
file2 = dir([p2fd,file2]);

load([p2fd,file1.name]);
load([p2fd,file2.name]);

%%
raw1 = save_data1{1};
raw2 = save_data2{1};

clear save_data*;

%%
cfg = [];

[avg1] = ft_timelockanalysis(cfg,raw1);%PAC
[avg2] = ft_timelockanalysis(cfg,raw2);%noPAC

%%
cfg                     = [];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.tapsmofrq           = 3*max(raw1.time{1});
cfg.output              = 'pow';
cfg.taper               = 'dpss';
cfg.keeptrials          = 'yes';

[freq1] = ft_freqanalysis(cfg,raw1);
[freq2] = ft_freqanalysis(cfg,raw2);

cfg = [];
cfg.frequency = [freq1.freq(1) 100];

[freq1] = ft_selectdata(cfg,freq1);

cfg = [];
cfg.frequency = [freq2.freq(1) 100];

[freq2] = ft_selectdata(cfg,freq2);

%%
x = mean(sqrt(abs(avg1.avg).^2),2);
z = (x-mean(x))./std(x);

idx = find(sign(z-2)==1);

hc = cell(length(idx),1);
for kt = 1:length(idx)
    hc{kt} = raw1.label{idx(kt)}(2);
end;

idxL = find(strcmp(hc(:),'L'));
idxR = find(strcmp(hc(:),'R'));

%%
figure;
subplot(221);
a = gca;
hold on;
h = [];
h(1) = plot(avg1.time,mean(avg1.avg(idxL,:),1).*1e15,'k','LineWidth',1);
h(2) = plot(avg1.time,mean(avg1.avg(idxR,:),1).*1e15,'r','LineWidth',1);
legend(h,'Left','Right');
legend('boxoff');
title('PAC');

subplot(222);
a = [a gca];
hold on;
plot(avg2.time,mean(avg2.avg(idxL,:),1).*1e15,'k');
plot(avg2.time,mean(avg2.avg(idxR,:),1).*1e15,'r');
title('no PAC');

axis(a,'tight');
lm1 = [min(min(avg1.avg([idxL;idxR],:))) max(max(avg1.avg([idxL;idxR],:)))];
lm2 = [min(min(avg2.avg([idxL;idxR],:))) max(max(avg2.avg([idxL;idxR],:)))];

lm = [min([lm1 lm2]) max([lm1 lm2])].*1e15;
set(a,'Xlim',[.15 .45]);
set(a,'XTick',[.2 .4]);
set(a(1),'YLim',[-max(max(lm)) max(max(lm))]./2);
set(a(2),'YLim',[-max(max(lm)) max(max(lm))]./2);
set(a(1),'YTick',round([-max(max(lm)) max(max(lm))]./4*10)/10);
set(a(2),'YTick',round([-max(max(lm)) max(max(lm))]./4*10)/10);
% lm = [-max(max(lm)) max(max(lm))];
% lm = lm.*10e12;
% lm = round(lm*10)/10;
% lm = lm*10e-12;
% set(a(1),'YTickLabel',lm);
% set(a(2),'YTickLabel',lm);

x = get(a(1),'YTick');
set(a,'YTick',[min(x) 0 max(x)]);
set(a,'FontName','Arial');
set(a,'Fontsize',12);
set(a,'LineWidth',3);
set(a,'ticklength',get(gca,'ticklength')*2);
set(gcf,'Color','w');

for it = 1:length(a)    
    xlabel(a(it),'Time [s]','Fontsize',14);
    ylabel(a(it),'Amplitude [fT]','Fontsize',14);
end;     

subplot(223);
hold on;
h = [];
h(1) = plot(freq1.freq,squeeze(mean(mean((freq1.powspctrm),2),1)).*1e15,'LineWidth',3);
h(2) = plot(freq2.freq,squeeze(mean(mean((freq2.powspctrm),2),1)).*1e15,'r','LineWidth',3);
axis tight;
ylabel('Power [fT^{2}/Hz]','Fontsize',14);
xlabel('Frequency [Hz]','Fontsize',14);
legend('PAC','no-PAC');
legend('boxoff');
% subplot(224);
% hold on;
% plot(sort(max(SNR1,[],1)));
% plot(sort(max(SNR2,[],1)),'r');
% axis tight;
% ylabel('Maximum SNR [dB]');
% xlabel('Trials');
x = get(gca,'YTick');
set(gca,'YTick',[0 max(x)]);
set(gca,'XTick',[10 50 90]);
set(a,'FontName','Arial');
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'ticklength',get(gca,'ticklength')*2);
set(gcf,'Color','w');

%%
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [8 12];
cfg.bpfilttype = 'fir';

[dum1] = ft_preprocessing(cfg,avg1);%PAC
[dum2] = ft_preprocessing(cfg,avg2);%noPAC

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [60 80];
cfg.bpfilttype = 'fir';

[dum3] = ft_preprocessing(cfg,avg1);%PAC
[dum4] = ft_preprocessing(cfg,avg2);%noPAC

cfg2  = [];
cfg2.method = 'template';
cfg2.layout = 'CTF275.lay';

cfg = [];
cfg.neighbours = ft_prepare_neighbours(cfg2);

dum1 = ft_megplanar(cfg,dum1);
dum1.avg = dum1.avg.^2;
dum1 = ft_combineplanar([],dum1);

dum2 = ft_megplanar(cfg,dum2);
dum2.avg = dum2.avg.^2;
dum2 = ft_combineplanar([],dum2);

dum3 = ft_megplanar(cfg,dum3);
dum3.avg = dum3.avg.^2;
dum3 = ft_combineplanar([],dum3);

dum4 = ft_megplanar(cfg,dum4);
dum4.avg = dum4.avg.^2;
dum4 = ft_combineplanar([],dum4);

cfg = [];
cfg.avgovertime = 'yes';

[dum1] = ft_selectdata(cfg,dum1);
[dum2] = ft_selectdata(cfg,dum2);
[dum3] = ft_selectdata(cfg,dum3);
[dum4] = ft_selectdata(cfg,dum4);

dum1.avg = sqrt(dum1.avg).*1e15;
dum2.avg = sqrt(dum2.avg).*1e15;
dum3.avg = sqrt(dum3.avg).*1e15;
dum4.avg = sqrt(dum4.avg).*1e15;

%%
cfg = [];
%cfg.grad = avg1.grad;
cfg.layout = 'CTF275.lay';
cfg.layout = ft_prepare_layout(cfg);
cfg.marker = 'off';
cfg.comment ='no';
%cfg.zlim = [0 1.5e-10];

figure;
subplot(221);
a = [a gca];
ft_topoplotER(cfg,dum3);%PAC
shading interp;lighting phong;
cb= colorbar;ca1 =caxis;
set(cb,'YTick',[min(get(cb,'YTick')) max(get(cb,'YTick'))]);
%set(cb,'YLim',round([ca1(1) ca1(2)].*100)/100);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
zlab = get(cb,'YLabel');
set(zlab,'String','RMS [fT]','Fontsize',14);
title('60-80Hz');
%set(cb,'YTick',round([ca1(1) ca1(2)].*10e13)/1000);
%set(cb,'YTickLabel',round([ca1(1) ca1(2)].*100)/100,'Fontsize',14);


subplot(222);
a = gca;
ft_topoplotER(cfg,dum4);%noPAC
shading interp;lighting phong;
cb = colorbar;caxis(ca1);
zlab = get(cb,'YLabel');
set(zlab,'String','RMS [fT]','Fontsize',14);
title('60-80Hz');
%set(cb,'YTick',round([ca(1) ca(2)].*100)/100);
%set(cb,'YTickLabel',round([ca(1) ca(2)].*100)/100,'Fontsize',14);
%set(cb,'YLim',round([ca(1) ca(2)].*100)/100);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
set(cb,'YTick',[min(get(cb,'YTick')) max(get(cb,'YTick'))]);

subplot(223);
a = [a gca];
ft_topoplotER(cfg,dum1);%PAC
shading interp;lighting phong;
cb = colorbar;ca2 = caxis;
zlab = get(cb,'YLabel');
set(zlab,'String','RMS [fT]','Fontsize',14);
title('8-12Hz');
% set(cb,'YTick',round([ca(1) ca(2)].*100)/100);
% set(cb,'YTickLabel',round([ca(1) ca(2)].*100)/100,'Fontsize',14);
% set(cb,'YLim',round([ca(1) ca(2)].*100)/100);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
set(cb,'YTick',[min(get(cb,'YTick')) max(get(cb,'YTick'))]);

subplot(224);
a = [a gca];
ft_topoplotER(cfg,dum2);%noPAC
title('8-12Hz');
shading interp;lighting phong;
cb = colorbar;caxis(ca2);
zlab = get(cb,'YLabel');
set(zlab,'String','RMS [fT]','Fontsize',14);
% set(cb,'YTick',round([ca(1) ca(2)].*100)/100);
% set(cb,'YTickLabel',round([ca(1) ca(2)].*100)/100,'Fontsize',14);
% set(cb,'YLim',round([ca(1) ca(2)].*100)/100);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
set(gcf,'Color','w');
set(cb,'YTick',[min(get(cb,'YTick')) max(get(cb,'YTick'))]);

%%
% savepath = '/home/rouxf/prj/TC/figures/';
% for it = 1:2
%     pos = get(gcf,'PaperPosition');
%     dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
%     set(gcf,'PaperPositionMode','auto');
%     if dim(2)>dim(1)
%         set(gcf,'PaperOrientation','landscape');
%     end;
%     it
%     print(it,'-dsvg','-opengl','-r600',[savepath,'MEGsensorLevel',num2str(it),'.svg']);
%     print(it,'-dpdf','-opengl','-r600',[savepath,'MEGsensorLevel',num2str(it),'.pdf']);
%     %plot2svg([savepath,'sourceModel_alpha_generator',num2str(it),'.svg'],it,'png');
% end;
return;
