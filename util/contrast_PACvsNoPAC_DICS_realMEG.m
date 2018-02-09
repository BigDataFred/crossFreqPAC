%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
%%
p2df = '~froux/froux/project_reaction_times/matFiles/';
files(1) = dir([p2df,'resting_state_recording_THA07_ICAcleanedMEGdata_eo.mat']);
files(2) = dir([p2df,'resting_state_recording_THA07_ICAcleanedMEGdata_ec.mat']);
%%
load([p2df,files(1).name]);
load([p2df,files(2).name]);

raw1 = dc1; clear dc1;
raw2 = dc2; clear dc2;

cfg = [];
cfg.channel = {'MEG'};

[dum1] = ft_selectdata(cfg,raw1);
[dum2] = ft_selectdata(cfg,raw2);

clear meg_data mri norm mni* loc r;
%%
grad1 = dum1.grad;
grad2 = dum2.grad;
%%
load([p2df,'individual_MNIwarped_grid.mat'],'grid','template_grid');
load([p2df,'individual_headmodel.mat'],'hdm');

Nx = length(template_grid.xgrid);
Ny = length(template_grid.ygrid);
Nz = length(template_grid.zgrid);
%%
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'powandcsd';
cfg.foilim = [10 10];
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 'maxperlen';

cfg.trials = 'all';
[csd1] = ft_freqanalysis(cfg,dum1);
[csd2] = ft_freqanalysis(cfg,dum2);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'powandcsd';
cfg.foilim = [70 70];
cfg.taper = 'dpss';
cfg.tapsmofrq = 10;
cfg.pad = 'maxperlen';

cfg.trials = 'all';
[csd3] = ft_freqanalysis(cfg,dum1);
[csd4] = ft_freqanalysis(cfg,dum2);

clear dum*;
%%
cfg         = [];
cfg.channel = {'MEG'};
cfg.grid    = grid;
cfg.headmodel     = hdm;
cfg.method  ='dics';
cfg.grid.dim=[Nx Ny Nz];
%cfg.normalize = 'yes';
cfg.keepleadfield = 'no';

cfg.dics.fixedori       = 'yes';
%cfg.dics.lambda = '1%'; %
cfg.dics.powmethod  = 'trace';
cfg.dics.projectnoise   = 'yes';
cfg.dics.keepfilter='yes';% these filters are for computing virtual electrodes later
cfg.dics.projectmom     = 'no';
cfg.dics.keepmom     = 'no';
cfg.dics.realfilter = 'yes';
cfg.dics.keepleadfield = 'no';

cfg.grad = grad1;
cfg.frequency = csd1.freq;
[source1] = ft_sourceanalysis(cfg,csd1);
source1 = rmfield(source1,'cfg');

cfg.grad = grad2;
cfg.frequency = csd2.freq;
[source2] = ft_sourceanalysis(cfg,csd2);
source2 = rmfield(source2,'cfg');

%align the beamformer with the grid positions
source1.pos = template_grid.pos;
source1.dim = template_grid.dim;

%align the beamformer with the grid positions
source2.pos = template_grid.pos;
source2.dim = template_grid.dim;

cfg         = [];
cfg.channel = {'MEG'};
cfg.grid    = grid;
cfg.headmodel     = hdm;
cfg.method  ='dics';
cfg.grid.dim=[Nx Ny Nz];
%cfg.normalize = 'yes';

cfg.dics.fixedori       = 'yes';
%cfg.dics.lambda = '1%'; %
cfg.dics.powmethod  = 'trace';
cfg.dics.projectnoise   = 'yes';
cfg.dics.keepfilter='yes';% these filters are for computing virtual electrodes later
cfg.dics.projectmom     = 'no';
cfg.dics.keepmom     = 'no';
cfg.dics.realfilter = 'yes';
cfg.dics.keepleadfield = 'no';

cfg.grad = grad1;
cfg.frequency = csd3.freq;
[source3] = ft_sourceanalysis(cfg,csd3);
source3 = rmfield(source3,'cfg');

cfg.grad = grad2;
cfg.frequency = csd4.freq;
[source4] = ft_sourceanalysis(cfg,csd4);
source4 = rmfield(source4,'cfg');

%align the beamformer with the grid positions
source3.pos = template_grid.pos;
source3.dim = template_grid.dim;

%align the beamformer with the grid positions
source4.pos = template_grid.pos;
source4.dim = template_grid.dim;
clear grad*;

%%
savepath = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';

savename1 = ['dics_PAC_sourceMap_alpha_eo_THA07.mat'];
savename2 = ['dics_noPAC_sourceMap_alpha_ec_THA07.mat'];
savename3 = ['dics_PAC_sourceMap_gamma_eo_THA07.mat'];
savename4 = ['dics_noPAC_sourceMap_gamma_ec_THA07.mat'];

save([savepath,savename1],'source1');
save([savepath,savename2],'source2');
save([savepath,savename3],'source3');
save([savepath,savename4],'source4');

%%
exit;