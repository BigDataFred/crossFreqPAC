%%
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
%%
p2d = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';
f1 = 'virtual_channels_1:150Hz_1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat';
f2 = '1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat';

load([p2d,f1]);
load([p2d,f2],'VCidx');

cfg                         = [];
cfg.channel                 = VC.label(VCidx);
cfg.pad                     = 'maxperlen';
cfg.method                  = 'mtmfft';
cfg.tapsmofrq               = 1;
cfg.taper                   = 'dpss';
cfg.output                  = 'fourier';

[phi]                       = ft_freqanalysis(cfg,VC);

[pow]                       = ft_freqdescriptives(cfg,phi);

cfg                         = [];
cfg.method                  = 'coh';
cfg.channel                 = VC.label(VCidx);
                            
[coh]                       = ft_connectivityanalysis(cfg,phi);
clear phi;
%% autoregressive model
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'no';

[dumds] = ft_resampledata(cfg,VC);

cfg = [];
cfg.channel = VC.label(VCidx);
cfg.order = 15;
cfg.toolbox = 'bsmart';

[mdata] = ft_mvaranalysis(cfg,dumds);
clear dumds;
%% transfer function
cfg = [];
cfg.method = 'mvar';

[mfreq] = ft_freqanalysis(cfg,mdata);
%%
cfg = [];
cfg.method = 'granger';
%cfg.granger.sfmethod = 'bivariate';

[granger] = ft_connectivityanalysis(cfg,mfreq);
%%
figure;
subplot(131);
a = gca;
plot(pow.freq,pow.powspctrm(1,:));%thalamus
subplot(132);
a = [a gca];
plot(pow.freq,pow.powspctrm(2,:));%FEF
subplot(133);
a = [a gca];
plot(pow.freq,pow.powspctrm(3,:));%left parietal cortex

for it = 1:length(a)
    axis(a(it),'tight');
    box(a(it),'off');
end;

set(a,'Xlim',[0 100]);
set(gcf,'Color','w');
set(a,'Fontsize',14);
%%
n1 = size(coh.cohspctrm,1);
n2 = size(coh.cohspctrm,2);
a =[];
figure;
c = 0;
for it = 1:n1
    for jt = 1:n2
        
        c = c+1;
        if it ~=jt
            subplot(n1,n2,c);
            %a(c) = gca;
            hold on;
            [a(c,:),h1,h2] = plotyy(coh.freq,squeeze(coh.cohspctrm(it,jt,:)),granger.freq,squeeze(granger.grangerspctrm(it,jt,:)));%coh
            set(h2,'Color','r');
        else
            subplot(n1,n2,c);
            a(c,:) = [gca gca];
            axis off;
        end;
        
    end;
end;

for it = 1:length(a)
    axis(a(it),'tight');
    box(a(it),'off');
end;
set(a(:,2),'YColor','r');
set(a,'Xlim',[0 100]);
set(a,'Ylim',[0 1]);
set(a,'YTick',[0 .5 1]);
set(gcf,'Color','w');
set(a,'Fontsize',14);
%%
cfg = [];
cfg.channel = VC.label(VCidx);% extract the gamma power for seed regions
cfg.bpfilter = 'yes';
cfg.bpfreq = [9 11];
cfg.bpfilttype = 'fir';
%cfg.hilbert = 'angle';% phase in rad
cfg.padding = max(VC.time{1});%
cfg.padtype = 'mirror';
cfg.continuous = 'no';

[alpha] = ft_preprocessing(cfg,VC);% single trial alpha activity

cfg = [];
cfg.channel = VC.label(VCidx);% extract the gamma power for seed regions
cfg.bpfilter = 'yes';
cfg.bpfreq = [65 75];
cfg.bpfilttype = 'fir';
%cfg.hilbert = 'abs'; % compute the power
cfg.padding = max(VC.time{1});%
cfg.padtype = 'mirror';
cfg.continuous = 'no';

[gamma] = ft_preprocessing(cfg,VC);
%%
n = VC.fsample*.1;
n = n*2+1;

x1 = zeros(length(VC.trial),n);
x2 = zeros(length(VC.trial),n);
x3 = zeros(length(VC.trial),n);
x4 = zeros(length(VC.trial),n);
x5 = zeros(length(VC.trial),n);
x6 = zeros(length(VC.trial),n);

for it = 1:length(VC.trial)
    
    y1 = alpha.trial{it}(1,:);%th
    y2 = alpha.trial{it}(2,:);%FEF
    y3 = alpha.trial{it}(3,:);%par
    
    [x1(it,:),lags] = xcorr(y1,y2,VC.fsample*.1,'coeff');
    [x2(it,:),~] = xcorr(y1,y3,VC.fsample*.1,'coeff');
    [x3(it,:),~] = xcorr(y2,y3,VC.fsample*.1,'coeff');
    
    y4 = gamma.trial{it}(1,:);%th
    y5 = gamma.trial{it}(2,:);%FEF
    y6 = gamma.trial{it}(3,:);%par
    
    [x4(it,:),~] = xcorr(y4,y5,VC.fsample*.1,'coeff');
    [x5(it,:),~] = xcorr(y4,y6,VC.fsample*.1,'coeff');
    [x6(it,:),~] = xcorr(y4,y5,VC.fsample*.1,'coeff');
end;