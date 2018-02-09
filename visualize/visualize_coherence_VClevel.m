%%
addpath('~froux/froux/fieldtrip-20151020');
ft_defaults;
%%
load 1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat VCidx;
%%
cfg= [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.channel =   VC.label(VCidx);
cfg.tapsmofrq = 1/max(VC.time{1});


[pow] = ft_freqanalysis( cfg , VC );

%%

cfg = [];
cfg.method = 'coh';

[coh] = ft_connectivityanalysis( cfg, pow );

%%

figure;
subplot(1,6,1);
a = gca;
plot(coh.freq,squeeze(coh.cohspctrm(1,2,:)));
subplot(1,6,2);
a = [a gca];
plot(coh.freq,squeeze(coh.cohspctrm(1,3,:)));


subplot(1,6,3);
a = [a gca];
plot(coh.freq,squeeze(coh.cohspctrm(2,1,:)));
subplot(1,6,4);
a = [a gca];
plot(coh.freq,squeeze(coh.cohspctrm(2,3,:)));

subplot(1,6,5);
a = [a gca];
plot(coh.freq,squeeze(coh.cohspctrm(3,1,:)));
subplot(1,6,6);
a = [a gca];
plot(coh.freq,squeeze(coh.cohspctrm(3,2,:)));

set(a,'XLim',[0 100]);
set(gcf,'Color','w');
for it = 1:length(a)
    xlabel(a(it),'Frequency (Hz)');
    ylabel(a(it),'Coherence (a.u.)');
end;