%%
p2d = '/home/rouxf/TC/matFiles/';
f1 = 'resting_state_recording_THA07_ICAcleanedMEGdataEC.mat';
f2 = 'resting_state_recording_THA07_GLMcleanedMEGdataEC.mat';

load([p2d,f1]);
load([p2d,f2]);

%%
cfg             = [];
cfg.derivative  = 'yes';

[dc2] = ft_preprocessing(cfg,dc2);

cfg             = [];
cfg.derivative  = 'yes';

[ec_data] = ft_preprocessing(cfg,ec_data);

%%
cfg             = [];
cfg.method      = 'mtmfft';
cfg.pad         = 'maxperlen';
cfg.output      ='pow';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 1/4;

[pow1] = ft_freqanalysis(cfg,dc2);
[pow2] = ft_freqanalysis(cfg,ec_data);

%%
figure;
hold on;
plot(pow1.freq,squeeze(mean(pow1.powspctrm,1)),'r');
plot(pow2.freq,squeeze(mean(pow2.powspctrm,1)));


