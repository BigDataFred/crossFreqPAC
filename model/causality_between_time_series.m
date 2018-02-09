%%
restoredefaultpath;
addpath('/home/rouxf/tbx/fieldtrip-20170618');
ft_defaults;
%%
Fs = 1.2e3;
t = 0:1/Fs:5;
%%
niter = 10;
v = zeros(niter,2);
idx1 = zeros(niter,2);
idx2 = zeros(niter,2);
granger = cell(niter,1);
c = zeros(niter,length(t)*2-1);
%%
for it = 1:niter
    
    orig = sin(2*pi*10.*t);
    x = orig;%randn(1,length(t));
    y = zeros(1,length(x));
    
    y = x;
    x = hilbert(x);
    x = x.*exp(1i.*pi/4);
    x = real(x);
    %%
%    for kt = 3:length(t)
 %       y(kt) = 0.9*y(kt-1)+0.5*x(kt-1)-0.8*y(kt-2);
 %   end;
    
    x = x+.1*randn(1,length(t));
    y = y+.1*randn(1,length(t));
    %%
    [sig1] = x;
    [sig2] = y;
    %%
    [c(it,:),lags] = xcorr(sig1,sig2,'coeff');
    %%
    dum = [];
    dum.label = {'dum_chan1','dum_chan2'};
    dum.trial{1} = [sig1;sig2];
    dum.time{1} = t;
    dum.fsample = Fs;
%     %%
%     w = 1;
%     cfg = [];
%     cfg.method = 'mtmconvol';
%     cfg.output = 'fourier';
%     cfg.taper = 'dpss';
%     cfg.tapsmofrq = 4;
%     cfg.pad = 'maxperlen';
%     cfg.padtype = 'zero';
%     cfg.toi = dum.time{1}(1)+(w/2+1e-2):0.025:dum.time{1}(end)-(w/2+1e-2);
%     cfg.foi = 1:2:100;
%     cfg.t_ftimwin = w*ones(1,length(cfg.foi));
%     
%     [freq] = ft_freqanalysis(cfg,dum);
    %% autoregressive model
    cfg = [];
    cfg.resamplefs = 200;
    cfg.detrend = 'no';
    
    [dumds] = ft_resampledata(cfg,dum);
    
    cfg = [];
    cfg.order = 15;
    cfg.toolbox = 'bsmart';
    
    [mdata] = ft_mvaranalysis(cfg,dumds);
    %% transfer function
    cfg = [];
    cfg.method = 'mvar';
    
    [mfreq] = ft_freqanalysis(cfg,mdata);
    %%
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.taper = 'dpss';
    cfg.tapsmofrq = .5;
    cfg.pad = 'maxperlen';
    cfg.padtype = 'zero';
    
    [pow] = ft_freqanalysis(cfg,dum);

    [~,idx1(it,1)] = max(pow.powspctrm(1,:));
    [~,idx1(it,2)] = max(pow.powspctrm(2,:));
    %%
    cfg = [];
    cfg.method = 'granger';
    %cfg.granger.sfmethod = 'bivariate';
    
    [granger{it}] = ft_connectivityanalysis(cfg,mfreq);
    
    
    [v(it,1),idx2(it,1)] = max(eval(['granger{',num2str(it),'}.',cfg.method,'spctrm(1,2,:)']));
    [v(it,2),idx2(it,2)] = max(eval(['granger{',num2str(it),'}.',cfg.method,'spctrm(2,1,:)']));

end;
%%
GC = zeros(niter,length(granger{1}.label),length(granger{1}.freq));
for it = 1:length(granger)
    GC(it,1,:) = squeeze(eval(['granger{',num2str(it),'}.',cfg.method,'spctrm(1,2,:)']))';
    GC(it,2,:) = squeeze(eval(['granger{',num2str(it),'}.',cfg.method,'spctrm(2,1,:)']))';
end;
GC = squeeze(mean(GC,1));
%%
figure;
subplot(321);
h = [];
hold on;
h(1) = plot(t,sig1);
h(2) = plot(t,sig2,'r');
xlim([0 .25]);ylim([-4.5 4.5]);
legend('S1','S2');
xlabel('Time [s]');
ylabel('Simulated amplitude [a.u.]');

subplot(322);
hold on;
h = [];
h(1) = plot(pow.freq,10*log10(pow.powspctrm(1,:)),'b-','LineWidth',2);
h(2) = plot(pow.freq,10*log10(pow.powspctrm(2,:)),'r');
plot([10 10],[min(min(10*log10(pow.powspctrm))) max(max(10*log10(pow.powspctrm)))],'k--');
axis tight;xlim([0 100]);
legend('S1','S2');
xlabel('Frequency [Hz]');
ylabel('Mean power [dB]');
set(gca,'XTick',[0:10:100]);

subplot(323);
hold on;
h = [];
h(1) = plot(granger{1}.freq,GC(1,:),'b');
h(2) = plot(granger{1}.freq,GC(2,:),'r');
plot([10 10],[min(min(GC)) max(max(GC))],'k--');
axis tight;
xlim([0 100]);
legend('S1->S2','S2->S1');
xlabel('Frequency [Hz]');
ylabel('Mean granger causality [a.u.]');

subplot(324);
hold on;
h = [];
h(1) = plot(granger{1}.freq(idx2(:,1)),v(:,1),'o','Color',[.25 .25 .25],'MarkerFaceColor','b');
h(2) = plot(granger{1}.freq(idx2(:,2)),v(:,2),'o','Color',[.25 .25 .25],'MarkerFaceColor','r');
plot([10 10],[min(min(v)) max(max(v))],'k--');
legend('S1->S2','S2->S1');
set(gca,'Color','w');
xlabel('Single trial granger peak frequency [Hz]');
ylabel('Single trial granger causality [a.u.]');
axis tight;
xlim([min(granger{1}.freq) 100]);
set(gca,'XTick',[0:10:100]);

subplot(325);
hold on;
plot(lags./Fs,mean(c,1),'k');box off;
d = lags(find(mean(c,1)==max(mean(c,1))))/Fs;
plot([d d],[min(mean(c,1)) max(mean(c,1))],'r--');
axis tight;
xlim([-.1 .1]);
ylabel('Correlation [a.u.]');
xlabel('Lag [ms]');

subplot(326);
l = zeros(size(c,1),1);
for it = 1:size(c,1)
    l(it) = lags(find(c(it,:) == max(c(it,:))))/Fs;
end;
[n,x] = hist(l,-.1:1/Fs:.1);
stem(x,n);box off;
axis tight;
ylabel('Count');
xlabel('Lag of max. correlation [ms]');

set(gcf,'Color','w');
%%
% [sig1] = sin(2*pi*10.*t);
% 
% [sig2] = sig1;
% [sig2] = hilbert(sig1);
% [sig2] = sig2.*exp(1i*(3/2*pi));
% [sig2] = real(sig2);
% 
% [sig1] = sig1+0.1*std(sig1).*randn(1,length(t));
% %[sig2] = sig1;
% %sig2 = fliplr(sig2);
% [sig2] = sig2+0.1*std(sig2).*randn(1,length(t));



