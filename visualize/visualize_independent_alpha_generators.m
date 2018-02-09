%%
[p2d] = '/home/rouxf/prj/TC/matFilesRev/';
[savepath] = '/home/rouxf/prj/TC/figures/';

load([p2d,'independent_alphaGenerators_timeSeries1.mat']);
load([p2d,'independent_alphaGenerators_coherenceMatrixANDspectra1.mat']);
load([p2d,'independent_alphaGenerators_amplitudeCorrelationMatrix1.mat']);
load([p2d,'independent_alphaGenerators_modulationIndexANDpah1.mat']);

%%
as = (.1/(2*pi));
X = params.phi1.*(180/pi);
Y = params.phi1.*as;

pbins = -pi:pi/8:pi;
phi = zeros(1,100);
for it = 1:100
    ix = randperm(length(params.phi1));
    ix = ix(1);
    phi(it) = params.phi1(ix);
end;
[n1,x1] = hist(phi,pbins);
[n2,x2] = hist(params.phi2,pbins);

figure;
subplot(121);
h1 = polar(pbins,n1);
set(gca,'LineWidth',3);
subplot(122);
h2 = polar(pbins,n2);
set(gca,'LineWidth',3);
lab1 = {'0', '90', '180','270','360'};
lab2 = {'\pi', '\pi/2 ','0', '-\pi2', '-\pi'};
lab3 = {'30', '60', '120','150','210','240','300','330'};

ix = findall(gcf,'type','text');
s = strtrim(get(ix,'String'));
for it = 1:length(lab1)
    set(ix(strcmp(s,lab1(it))),'String',lab2(it),'Fontsize',14);
end;
for it = 1:length(lab3)
    set(ix(strcmp(s,lab3(it))),'String','');
end;
set(h1,'LineWidth',3,'Color','r');
set(h2,'LineWidth',3,'Color','b');
set(gca,'LineWidth',3);
set(gcf,'Color','w');
set(gca,'FontName','Arial')

%%
figure;
subplot(331);
b=gca;
hold on;
h =[];
h(1) = plot(params.t,squeeze(mean(dat(1,:,:),2)));
h(2) =plot(params.t,.06+squeeze(mean(dat(3,:,:),2)),'r');
sc1 = mean(squeeze(mean(dat(1,:,:),2)));
sc2 = mean(.06+squeeze(mean(dat(3,:,:),2)));
set(gca,'YTick',[sc1-.01 sc1+.01 sc2-.01 sc2+.01]);
set(gca,'YTickLabel',[-0.01 0.01 -0.01 0.01]);
yls = get(gca,'YLim');
legend(h,'S1','S3');
legend('boxoff');
axis tight;
subplot(332);
a = [gca];
hold on;
h= [];
h(1) = plot(f,squeeze(S1(:,1,1)),'LineWidth',2);
title('S1');
subplot(333);
a = [a gca];
hold on;
h = [];
h(1) = plot(f,squeeze(S1(:,3,3)),'r','LineWidth',2);
title('S3');

subplot(334);
b=[b gca];
hold on;
h =[];
h(1) = plot(params.t,squeeze(mean(dat(2,:,:),2)),'Color',[.75 .75 .75]);
h(2) =plot(params.t,.06+squeeze(mean(dat(4,:,:),2)),'Color',[0 0 0]);
sc1 = mean(squeeze(mean(dat(2,:,:),2)));
sc2 = mean(.06+squeeze(mean(dat(4,:,:),2)));
set(gca,'YTick',[sc1-.01 sc1+.01 sc2-.01 sc2+.01]);
set(gca,'YTickLabel',[-0.01 0.01 -0.01 0.01]);
yls = [yls get(gca,'YLim')];
legend(h,'S2','S4');
legend('boxoff');
axis tight;
subplot(335);
a = [a gca];
hold on;
h = [];
h(1) = plot(f,squeeze(S1(:,2,2)),'Color',[.75 .75 .75],'LineWidth',2);
title('S2');
subplot(336);
a = [a gca];
hold on;
h = [];
h(1) = plot(f,squeeze(S1(:,4,4)),'k','LineWidth',2);
title('S4');

subplot(337);
b=[b gca];
hold on;
h =[];
h(1) = plot(params.t,squeeze(mean(dat(5,:,:),2)),'Color',[0 .5 .5]);
h(2) =plot(params.t,.3+squeeze(mean(dat(6,:,:),2)),'Color',[.9 .5 .5]);
sc1 = mean(squeeze(mean(dat(5,:,:),2)));
sc2 = mean(.3+squeeze(mean(dat(6,:,:),2)));
set(gca,'YTick',[sc1-.1 sc1+.1 sc2-.1 sc2+.1]);
set(gca,'YTickLabel',[-0.1 0.1 -0.1 0.1]);
yls = [yls get(gca,'YLim')];
legend(h,'S5','S6');
legend('boxoff');
axis tight;
subplot(338);
a = [a gca];
hold on;
h = [];
h(1) = plot(f,squeeze(S1(:,5,5)),'Color',[0 .5 .5],'LineWidth',2);
title('S5');
subplot(339);
a = [a gca];
hold on;
h = [];
h(1) = plot(f,squeeze(S1(:,6,6)),'Color',[.9 .5 .5],'LineWidth',2);
title('S6');

yl = [];
for it  =1:length(a)
    xlabel(a(it),'Frequency [Hz]','Fontsize',14);
    ylabel(a(it),'Power [a.u.]','Fontsize',14);
    yl(it,:) = get(a(it),'YLim');
end;
%set(a,'YLim',[min(min(yl)) max(max(yl))]);
set(a,'XTick',[10 50 90]);
set(a,'FontName','Arial');

set(b,'XLim',[.15 .45]);
set(b,'XTick',[.2 .4]);
set(b(1:2),'YLim',[min(min(yls(1:4))) max(max(yls(1:4)))]);
%set(b(3),'YLim',[min(min(yls(5:6))) max(max(yls(5:6)))]);
set(b(3),'YLim',[-0.15 max(max(yls(5:6)))]);
%set(b(2),'YTick',[-0.01 0 0.01]);
%set(b([1 3]),'YTick',[-0.1 0 0.1]);
for it  =1:length(b)
    xlabel(b(it),'Time [s]','Fontsize',14);
    ylabel(b(it),'Amplitude [nA*m]','Fontsize',14);    
end;
set(gcf,'Color','w');
set(a,'Fontsize',12);
set(a,'LineWidth',3);
set(a,'ticklength',get(gca,'ticklength')*2);
set(b,'Fontsize',12);
set(b,'LineWidth',3);
set(b,'ticklength',get(gca,'ticklength')*2);
set(b,'FontName','Arial');

%%
figure;
imagesc(r);
axis xy;
caxis([-1 1]);
set(gca,'XTick',1:6);
set(gca,'YTick',1:6);
set(gca,'XTickLabel',{'S1' 'S2' 'S3' 'S4' 'S5' 'S6'},'Fontsize',14);
set(gca,'YTickLabel',{'S1' 'S2' 'S3' 'S4' 'S5' 'S6'},'Fontsize',14);
title('Amplitude correlation matrix');
set(gcf,'Color','w');
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'ticklength',get(gca,'ticklength')*2);
set(gca,'FontName','Arial');

figure;
colorbar;
caxis([-1 1]);
axis off;
cb= colorbar;
zlab = get(cb,'Ylabel');
set(zlab,'String','Spearman Coefficient [a.u.]','Fontsize',14);
set(cb,'YTick',[-1 1]);
set(gcf,'Color','w');
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'ticklength',get(gca,'ticklength')*2);
set(cb,'FontName','Arial');

%%
idx = 1:size(C,2);
figure;
cnt = 0;
a = zeros(1,size(C,2)*size(C,3));
for it = 1:size(C,2)
    for jt = 1:length(idx)
        
        if it ==idx(jt)
        else
            cnt = cnt+1;
            subplot(5,5,cnt);
            a(cnt) = gca;
            plot(f,squeeze(C(:,it,idx(jt))),'LineWidth',3);
            axis tight;
            ylim([0 1]);
            title(['S',num2str(it),'-S',num2str(idx(jt))]);
        end;
    end;
    idx(1) = [];
    cnt = cnt+it;
end;
a(a==0) = [];

set(a,'XLim',[0 100]);
set(a,'YLim',[0 1.05]);
set(a,'YTick',[0 1]);
set(a,'XTick',[10 50 90]);
set(a,'Box','on');
for it = 6%[1 6 11 16]
    xlabel(a(it),'Frequency [Hz]','Fontsize',14);
    ylabel(a(it),'Coherence [a.u.]','Fontsize',14);
end;

set(gcf,'Color','w');
set(a,'Fontsize',12);
set(a,'LineWidth',3);
set(a,'ticklength',get(gca,'ticklength')*2);
set(a,'FontName','Arial');

%%
pbins2 = pbins.*(180/pi)+180;

selIx1 = find(pfoi >=9 & pfoi <=11);
selIx2 = find(afoi >=60 & afoi <=90);

x = [];
for it = 1:length(PAH)
    x(it,:) = squeeze(mean(mean(mean(PAH{it}(:,selIx1,selIx2,:),3),2),1));    
end;

figure;
cnt = 0;
cnt2 = size(MI,1);
b = zeros(1,size(MI,1));
a = zeros(1,size(MI,1));
for it = 1:size(MI,1)
    
    cnt = cnt+1;
    subplot(2,6,cnt);
    b(it) = gca;
    imagesc(pfoi,afoi,squeeze(mean(MI(it,:,:,:),2))');
    if it ==5
        ca = caxis;
    end;
    axis xy;
    title(['S',num2str(it)]);
    
    cnt2 = cnt2+1;
    subplot(2,6,cnt2);
    a(it) =gca;
    hold on;
    if ismember(it,[3 4])
        y = cos(pbins)*max(x(it,:));
    elseif ismember(it,[5 6]);
        y = sin(pbins)*max(x(it,:));
    else
        y = mean(x(it,:))*ones(1,length(pbins));      
    end;    
    bar([pbins2 pbins2+361],[x(it,:) x(it,:)]);
    %plot([pbins2 pbins2+361],[y y],'r','LineWidth',2);
    axis tight;yl(1,:) = get(gca,'YLim');

end;

for it = 1:length(b)
    %caxis(b(it),ca);
end;

for it = 1:length(b)
    xlabel(b(it),'Frequency [Hz]','Fontsize',14);
    ylabel(b(it),'Frequency [Hz]','Fontsize',14);
end;
set(b,'YTick',[afoi(1):40:afoi(end)]);
set(b,'XTick',[pfoi(1):8:pfoi(end)]);

for it = 1:length(a)
    xlabel(a(it),'\alpha-phase [rad]','Fontsize',14);
    ylabel(a(it),'\gamma-power [a.u.]','Fontsize',14);
end;
set(a,'XTick',[0 180 360 540 720]);
%set(a,'XTickLabel',{'-\pi' '-\pi/2' '0' '\pi/2' '\pi'},'Fontsize',14);
%set(a,'YTick',[0 1]);
%set(a,'YTickLabel',[0 1],'Fontsize',14);
%set(a,'YLim',[min(min(yl)) max(max(yl))]);
set(a,'Fontsize',12);
set(a,'LineWidth',3);
set(a,'ticklength',get(gca,'ticklength')*2);
set(a,'FontName','Arial');

set(b,'FontName','Arial');
set(b,'Fontsize',12);
set(b,'LineWidth',3);
set(b,'ticklength',get(gca,'ticklength')*2);
set(gcf,'Color','w');

%%
figure;
colorbar;
caxis(ca);
axis off;
cb= colorbar;
zlab = get(cb,'Ylabel');
set(zlab,'String','Modulation Index [10^{-3}]','Fontsize',14);
set(cb,'YTick',round([ca(1) ca(2)].*1e4)/1e4);
set(cb,'YTickLabel',(round([ca(1) ca(2)].*1e4)/1e4).*1e3,'Fontsize',14);
set(cb,'YLim',round([ca(1) ca(2)].*1e4)/1e4);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',max(get(gca,'ticklength'))*2);
set(gcf,'Color','w');

% %%
% cnt = 4;
% for it =2%[1 2 3 4 5 7 9]
%     cnt = cnt+1;
%     pos = get(gcf,'PaperPosition');
%     dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
%     set(gcf,'PaperPositionMode','auto');
% %     if dim(2)>dim(1)
% %         set(gcf,'PaperOrientation','landscape');
% %     end;
%     it
%     print(it,'-dsvg','-opengl','-r800',[savepath,'sourceModel_alphaGenerator',num2str(it),'.svg']);
%     print(it,'-dpdf','-opengl','-r800',[savepath,'sourceModel_alphaGenerator',num2str(it),'.pdf']);
%     %plot2svg([savepath,'sourceModel_alpha_generator',num2str(it),'.svg'],it,'png');
% end;
