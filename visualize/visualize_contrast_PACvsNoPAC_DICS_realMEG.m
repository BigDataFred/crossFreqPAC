%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
p2df = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';

%% normalize individual T1 to template brain
load([p2df,'THA07_V2.mri_MRI4MNI.mat']);

cfg = [];
cfg.spmversion  = 'spm8';
cfg.template = '/bcbl/home/home_a-f/froux/spm8 2/templates/T1.nii';
cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);


%%
load([p2df,'dics_PAC_sourceMap_alpha_eo_THA07.mat']);
load([p2df,'dics_noPAC_sourceMap_alpha_ec_THA07.mat']);
load([p2df,'dics_PAC_sourceMap_gamma_eo_THA07.mat']);
load([p2df,'dics_noPAC_sourceMap_gamma_ec_THA07.mat']);

dum1 = source1;
dum1.avg.pow = 20*log10(source2.avg.pow./source1.avg.pow);

dum2 = source3;
dum2.avg.pow = 20*log10(source4.avg.pow./source3.avg.pow);

cfg = [];
cfg.parameter = 'all';

[int1] = ft_sourceinterpolate(cfg,dum1,norm);
[int2] = ft_sourceinterpolate(cfg,dum2,norm);
%%
cfg = [];
cfg.method          = 'slice';
cfg.funparameter    = 'pow';
cfg.maskparameter   = cfg.funparameter;
cfg.opacitymap      = 'rampup';
cfg.funcolormap     = 'jet';
cfg.colorbar        = 'no';
cfg.nslices         = 40;

cfg.opacitylim = [max(dum1.avg.pow)*0.001 max(dum1.avg.pow)];
ft_sourceplot(cfg,int1);

%%
cfg.opacitymap      = 'rampdown';
cfg.opacitylim = [min(dum2.avg.pow) min(dum2.avg.pow)*0.001];
ft_sourceplot(cfg,int2);




