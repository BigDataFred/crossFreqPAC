restoredefaultpath;
addpath('/media/rouxf/rds-share/Common/fieldtrip-20170115/');           
ft_defaults;
savepath = '~rouxf/TC/matFiles/'; 

%%
cfg = [];
cfg.dataset     = '~rouxf/TC/THA07_DEVELOPMENT_20090924_12.ds';
cfg.continous   = 'yes';

cfg.bpfilter    = 'yes';
cfg.bpfreq      = [1 200];
cfg.bpfilttype  = 'fir';

cfg.dftfilter   = 'yes';
cfg.dftfreq     = [50:50:200];

cfg.padding     = 20;
cfg.padtype     = 'data';

cfg.channel = {'MEG'};

cfg.trl         = [20*1200 220*1200 0];
[eo_data] = ft_preprocessing(cfg);

cfg.trl         = [260*1200 460*1200 0];
[ec_data] = ft_preprocessing(cfg);

%%
cfg             = [];
cfg.length      = [ 4 ];

[eo_data]      = ft_redefinetrial(cfg,eo_data);
[ec_data]      = ft_redefinetrial(cfg,ec_data);

%%
cfg             = [];
cfg.resamplefs  = 150;
cfg.detrend     = 'no';

[ds1] = ft_resampledata(cfg,eo_data);
[ds2] = ft_resampledata(cfg,ec_data);

%%
cfg             = []; 
cfg.method      ='runica';

[comp1] = ft_componentanalysis(cfg,ds1);
[comp2] = ft_componentanalysis(cfg,ds2);

%%
cfg             = [];
cfg.method      ='runica';
cfg.unmixing    = comp1.unmixing;

[comp1] = ft_componentanalysis(cfg, eo_data);

cfg             = [];
cfg.method      ='runica';
cfg.unmixing    = comp2.unmixing;

[comp2] = ft_componentanalysis(cfg, ec_data);

%%
cfg = [];
cfg.dataset     = '~rouxf/TC/THA07_DEVELOPMENT_20090924_12.ds';
cfg.continous   = 'yes';

cfg.bpfilter    = 'yes';
cfg.bpfreq      = [1 200];
cfg.bpfilttype  = 'fir';

cfg.dftfilter   = 'yes';
cfg.dftfreq     = [50:50:200];

cfg.padding     = 20;
cfg.padtype     = 'data';

cfg.channel = {'EEG057','EEG058'};

cfg.trl         = [20*1200 220*1200 0];
[eeg_data1] = ft_preprocessing(cfg);

cfg.trl         = [260*1200 460*1200 0];
[eeg_data2] = ft_preprocessing(cfg);

%%
cfg             = [];
cfg.length      = [ 4 ];

[eeg_data1]      = ft_redefinetrial(cfg,eeg_data1);
[eeg_data2]      = ft_redefinetrial(cfg,eeg_data2);

%%
addpath('~rouxf/icacleaning/');

[ecg_idx1,qrs1] = clean_ecg(comp1);
[eog_idx1,eog1,~] = clean_eog(comp1,eeg_data1);

[ecg_idx2,qrs2] = clean_ecg(comp2);
[eog_idx2,eog2,~] = clean_eog(comp2,eeg_data2);

%%
% cfg = [];
% cfg.layout = 'CTF275.lay';
% cfg.channel = [eog_idx1 ecg_idx1];
% 
% ft_databrowser(cfg,comp1);

%%
cfg                 = [];
cfg.component       = [ecg_idx1 eog_idx1];

[dc1] = ft_rejectcomponent(cfg,comp1,eo_data);

cfg                 = [];
cfg.component       = [ecg_idx2 eog_idx2];

[dc2] = ft_rejectcomponent(cfg,comp2,ec_data);

%%
[eo_data] = glmcleaning(eo_data,ecg_idx1,qrs1);
out_name = 'resting_state_recording_THA07_GLMcleanedMEGdataEO.mat';
save([savepath,out_name],'eo_data');

[ec_data] = glmcleaning(ec_data,ecg_idx2,qrs2);
out_name = 'resting_state_recording_THA07_GLMcleanedMEGdataEC.mat';
save([savepath,out_name],'ec_data');

%%
cfg                 = [];
cfg.layout          = 'CTF275.lay';

[lay] = ft_prepare_layout(cfg);

cfg2                = [];
cfg2.layout         ='CTF275.lay';
cfg2.method         ='template';
cfg2.grad           = dc1.grad;

cfg                 = [];
cfg.method          = 'average';
cfg.missingchannel  = setdiff(lay.label(1:275),dc1.label);
cfg.neighbours      = ft_prepare_neighbours(cfg2); 

[dc1] = ft_channelrepair(cfg,dc1);
out_name = 'resting_state_recording_THA07_ICAcleanedMEGdataEO.mat';
save([savepath,out_name],'dc1');

%%
cfg                 = [];
cfg.layout          = 'CTF275.lay';

[lay] = ft_prepare_layout(cfg);

cfg2                = [];
cfg2.layout         ='CTF275.lay';
cfg2.method         ='template';
cfg2.grad           = dc2.grad;

cfg                 = [];
cfg.method          = 'average';
cfg.missingchannel  = setdiff(lay.label(1:275),dc2.label);
cfg.neighbours      = ft_prepare_neighbours(cfg2); 

[dc2] = ft_channelrepair(cfg,dc2);
out_name = 'resting_state_recording_THA07_ICAcleanedMEGdataEC.mat';
save([savepath,out_name],'dc2');



