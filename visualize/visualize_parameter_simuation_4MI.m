%%
tl = 1:30;

Fs = 300:100:1200;
pfoi = 4:2:20;
afoi = 45:5:140;

pbins = -pi:pi/8:pi;

fn = 'parameters4MI_alphaGammaCoupling.mat';
p2df  = '/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/matFiles/';
load([p2df,fn]);
%%
tl = 1:30;

Fs = 300:100:1200;
pfoi = 4:2:20;
afoi = 45:5:140;

pbins = -pi:pi/8:pi;
%%
aix = find(afoi  >= 65 & afoi <=75);
pix = find(pfoi  >= 9 & pfoi <=11);

figure;
subplot(2,3,[1]);
hold on;
for it = 1:length(tl)
    X = Fs;
    Y = squeeze(mean(mean(MI(:,it,pix,aix),4),3));
    if it == 1
        %plot(X,Y,'s-','Color',[.75 .75 .75],'MarkerFaceColor',[.75 .75 .75]);
    elseif it ==length(tl)
        %plot(X,Y,'o-','Color',[.75 .75 .75],'MarkerFaceColor',[.75 .75 .75]);
    else
        %plot(X,Y,'Color',[.75 .75 .75]);
    end;
end;
X = Fs;
Y = squeeze(mean(mean(mean(MI(:,:,pix,aix),4),3),2));
Y2 = squeeze(std(mean(mean(MI(:,:,pix,aix),4),3),0,2));

plot(X,Y,'Color',[.9 0 0],'LineWidth',3);
plot(X,Y-Y2,'--','Color',[.75 .75 .75],'LineWidth',3);
plot(X,Y+Y2,'--','Color',[.75 .75 .75],'LineWidth',3);

axis tight;
ylabel('Modulation index [a.u.]');
xlabel('Sampling frequency [Hz]');

subplot(2,3,3);
a = gca;
Z = squeeze(MI(1,1,:,:));

pcolor(pfoi,afoi,Z');
ca = caxis;
shading interp;
lighting phong;
cb = colorbar;

subplot(2,3,6);
a = [a gca];
Z = squeeze(MI(end,end,:,:));
pcolor(pfoi,afoi,Z');
shading interp;
lighting phong;
caxis(ca);
cb = colorbar;

set(gcf,'Color','w');

for it = 1:length(a)
    xlabel(a(it),'Frequency for phase [Hz]');
    ylabel(a(it),'Frequency for amplitude [Hz]');
end;
%%
[MI] = zeros(length(pfoi),length(afoi));

[t] = 0:1/Fs(end):tl(end);

[asig] = sin(2*pi*10.*t);%
[gsig] = (sin(2*pi*10.*t)+1).*sin(2*pi*70.*t);

[sig] = asig+gsig+1;%
[noise] = (mean(sig)+2*std(sig).*randn(1,length(t)));

[sig] = sig + noise;
%%
dum = [];
dum.label = {'dummy_channel'};
dum.cfg = [];
dum.trial{1} = sig;
dum.time{1} = t;

if matlabpool('size') ==0
    matlabpool(length(pfoi));
end;

parfor kt = 1:length(pfoi)
    
    cfg = [];
    cfg.hilbert = 'angle';
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [pfoi(kt)-1 pfoi(kt)+1];
    cfg.bpfilttype = 'fir';
    
    [phi] = ft_preprocessing(cfg,dum);
    p = phi.trial{1};
    
    [mi] = zeros(length(afoi),1);
    [X]  = zeros(length(afoi),length(pbins));

    for mt = 1:length(afoi)
        
        cfg = [];
        cfg.hilbert = 'abs';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [afoi(mt)-5 afoi(mt)+5];
        cfg.bpfilttype = 'fir';
        
        [amp] = ft_preprocessing(cfg,dum);
        a = amp.trial{1}.^2;
        
        for zt = 1:length(pbins)-1
            a2 = a;
            pbins2 = pbins;
            pidx = find( p >= pbins2(zt) & p < pbins2(zt+1));
            
            X(mt,zt) = mean(a2(pidx));
            
        end;
        X(mt,end) = X(mt,1);
        
        n = length(pbins);
        X(mt,:) = X(mt,:)./sum(X(mt,:),2);
        H = -sum(log(X(mt,:)).*X(mt,:),2);
        mi(mt) = (log(n)-H)./log(n);
        
    end;
    PAH(kt,:,:) = X;
    MI(kt,:) = mi;
end;
matlabpool close;
%%
subplot(2,3,4);
plot(pbins,squeeze(mean(PAH(pix,aix,:),2)),'bo-','LineWidth',3);
axis tight;
xlabel('Alpha phase [rad]');
ylabel('Gamma power [a.u.]');